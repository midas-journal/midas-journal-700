/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDicomDirectoryToImage.txx,v $
  Language:  C++
  Date:      $Date: 2009/11/01 21:06:52 $
  Version:   $Revision: 1.1 $
=========================================================================*/
#ifndef _itkDicomDirectoryToImage_txx
#define _itkDicomDirectoryToImage_txx
#include "itkDicomDirectoryToImage.h"

// using gdcm 2.x
#include "gdcmAttribute.h"
#include "gdcmDataSet.h"
#include "gdcmImage.h"
#include "gdcmImageReader.h"
#include "gdcmReader.h"
#include "gdcmWriter.h"
#include "itkDirectory.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVersion.h"
#include "itksys/SystemTools.hxx"

#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <map>

#include <sys/stat.h>

namespace itk
{

// constructor
template<class PT>
DicomDirectoryToImage<PT>
::DicomDirectoryToImage(const std::string dirname=""):
    m_Verbose(false),
    m_DicomDirectory(dirname),
    m_OutputDirectory(""),
    m_MatchOnSeries(""),
    m_InterpretedCriteria(""),
    m_MinimumSlices(0),
    m_Anonymize(false),
    m_CSVMappingFile(""),
    m_UseAccessionNumberInNaming(false),
    m_FileNameExtension("tif"),
    m_CriteriumSingle(""),
    m_Recursive(false)
{
    m_CriteriumAnd.clear();
    m_CriteriumOr.clear();
    m_DirsToImport.clear();
    m_LogMessages.clear();
}
//destructor
template<class PT>
DicomDirectoryToImage<PT>
::~DicomDirectoryToImage()
{
}
 
// interpret match criterium
template<class PT>
bool DicomDirectoryToImage<PT>
::InterpretMatchCriterium() 
{

    // first check if string is given
    if (m_MatchOnSeries.empty())
    {
        return false;
    }

    // check whether AND/OR operators are used. If so, split
    // up the string otherwise a single matchcriterium is used
    if ( m_MatchOnSeries.find(" AND ") == std::string::npos &&
         m_MatchOnSeries.find(" OR ") == std::string::npos)
    {
        // one string only
        m_CriteriumSingle = m_MatchOnSeries;
        std::transform (m_CriteriumSingle.begin(),
                        m_CriteriumSingle.end(),
                        m_CriteriumSingle.begin(),
                        tolower);
    }
    else
    {
        // more than one string
        // split this string with space as delimiter
        std::string                     buffer;
        std::vector <std::string>       matchcriteria;
        std::stringstream               ss(m_MatchOnSeries);

        while (ss >> buffer)
            matchcriteria.push_back(buffer);

        // first and last element should not equal AND/OR
        if (matchcriteria[0] == "AND" ||
            matchcriteria[0] == "OR" ||
            matchcriteria[matchcriteria.size()-1] == "AND" ||
            matchcriteria[matchcriteria.size()-1] == "OR")
        {
            if (m_Verbose)
            {
                std::cout << "invalid series match criterium" << std::endl;
            }
            return false;
        }

        // traverse all strings
        std::vector<bool> visited(matchcriteria.size());

        // obtain arguments of all ANDs
        for (unsigned int i=1; i<matchcriteria.size()-1; ++i)
        {
            if (matchcriteria[i] == "AND" && !visited[i])
            {
                visited[i] = 1;

                // process left from AND
                unsigned int j = i-1;
                if (matchcriteria[j] == "AND" || matchcriteria[j] == "OR")
                {
                    if (m_Verbose)
                    {
                        std::cout << "invalid series match criterium" << std::endl;
                    }
                    return false;
                }
                std::string left("");
                bool finished = false;
                while (!finished)
                {
                    if ( !visited[j] &&
                         matchcriteria[j] != "AND" &&
                         matchcriteria[j] != "OR" )
                    {
                        visited[j] = 1;
                        if (left.empty())
                            left = matchcriteria[j];
                        else
                            left = matchcriteria[j] + " " + left;

                        if (j==0)
                            finished = true;
                        else 
                            --j;
                    }
                    else finished = true;
                }
                if (!left.empty())
                {
                    std::transform (left.begin(), left.end(), left.begin(), tolower);
                    m_CriteriumAnd.push_back(left);
                }

                // process right from AND
                j = i + 1;
                if (matchcriteria[j] == "AND" || matchcriteria[j] == "OR" )
                {
                    if (m_Verbose)
                    {
                        std::cout << "invalid series match criterium" << std::endl;
                    }
                    return false;
                }

                finished = false;
                std::string right ("");
                while (!finished)
                {
                    if ( !visited[j] &&
                         matchcriteria[j] != "AND" &&
                         matchcriteria[j] != "OR" )

                    {
                        visited[j] = 1;

                        if (right.empty())
                            right = matchcriteria[j];
                        else
                            right = right + " " + matchcriteria[j];

                        if ( j == matchcriteria.size() - 1)
                            finished = true;
                        else 
                            ++j;
                    }
                    else finished = true;
                }
                if (!right.empty())
                {
                    std::transform (right.begin(), right.end(), right.begin(), tolower);
                    m_CriteriumAnd.push_back(right);
                }
            }
        } // end arguments of ANDs

        // obtain arguments of all ORs 
        for (unsigned int i=1; i<matchcriteria.size()-1; ++i)
        {
            if (matchcriteria[i] == "OR" && !visited[i])
            {
                visited[i] = 1;

                // process left from OR
                unsigned int j = i-1;
                if (matchcriteria[j] == "AND" || matchcriteria[j] == "OR")
                {
                    if (m_Verbose)
                    {
                        std::cout << "invalid series match criterium" << std::endl;
                    }
                    return false;
                }
                std::string left("");
                bool finished = false;
                while (!finished)
                {
                    if ( !visited[j] &&
                         matchcriteria[j] != "AND" &&
                         matchcriteria[j] != "OR" )
                    {
                        visited[j] = 1;
                        if (left.empty())
                            left = matchcriteria[j];
                        else
                            left = matchcriteria[j] + " " + left;

                        if (j==0)
                            finished = true;
                        else 
                            --j;
                    }
                    else finished = true;
                }
                if (!left.empty())
                {
                    std::transform (left.begin(), left.end(), left.begin(), tolower);
                    m_CriteriumOr.push_back(left);
                }

                // process right from AND
                j = i + 1;
                if (matchcriteria[j] == "AND" || matchcriteria[j] == "OR" )
                {
                    if (m_Verbose)
                    {
                        std::cout << "invalid series match criterium" << std::endl;
                    }
                    return false;
                }

                finished = false;
                std::string right ("");
                while (!finished)
                {
                    if ( !visited[j] &&
                         matchcriteria[j] != "AND" &&
                         matchcriteria[j] != "OR" )

                    {
                        visited[j] = 1;

                        if (right.empty())
                            right = matchcriteria[j];
                        else
                            right = right + " " + matchcriteria[j];

                        if ( j == matchcriteria.size() - 1)
                            finished = true;
                        else 
                            ++j;
                    }
                    else finished = true;
                }
                if (!right.empty())
                {
                    std::transform (right.begin(), right.end(), right.begin(), tolower);
                    m_CriteriumOr.push_back(right);
                }
            }
        } // end arguments of Or

    } // end else multiple strings

    return true;
}
// print match criterium
template<class PT>
void DicomDirectoryToImage<PT>
::PrintMatchCriterium() 
{
    if (!m_CriteriumSingle.empty())
    {
        m_InterpretedCriteria = m_CriteriumSingle;
        if (m_Verbose)
        {
            std::cout << " series criterium    \t" << m_CriteriumSingle << std::endl;
        }
    }

    if (!m_CriteriumAnd.empty())
    {
        m_InterpretedCriteria = m_CriteriumAnd[0] + " ";

        if (m_Verbose)
        {
            std::cout << " series criteria     \t(" << m_CriteriumAnd[0] << " ";
        }

        for (unsigned int i=1; i<m_CriteriumAnd.size(); ++i)
        {
            m_InterpretedCriteria = m_InterpretedCriteria + "AND " +
                                    m_CriteriumAnd[i];
            if (m_Verbose)
                std::cout << "AND " << m_CriteriumAnd[i];
        }
        if (m_Verbose)
        {
            std::cout << ")";
        }
    }

    if (!m_CriteriumOr.empty())
    {
        if (m_CriteriumAnd.empty())
        {
            if (m_Verbose)
            {
                std::cout << " series criteria     \t";
            }
        }
        else
        {
            m_InterpretedCriteria = m_InterpretedCriteria + " OR ";
            if (m_Verbose)
            {
                std::cout << " OR (";
            }
        }

        m_InterpretedCriteria = m_InterpretedCriteria + m_CriteriumOr[0] + " ";

        if (m_Verbose)
        {
            std::cout << m_CriteriumOr[0] << " ";
        }

        for (unsigned int i=1; i<m_CriteriumOr.size(); ++i)
        {
            m_InterpretedCriteria = m_InterpretedCriteria + "OR " +
                                    m_CriteriumOr[i];
            if (m_Verbose)
                std::cout << "OR " << m_CriteriumOr[i];
        }
        if (m_Verbose)
        {
            std::cout << ")";
        }
            

    }
    if (!m_CriteriumAnd.empty() || !m_CriteriumOr.empty())
    {
        if (m_Verbose)
        {
            std::cout << std::endl;
        }
    }
}
// getalldicomseriesfromsubdir
template<class PT>
typename DicomDirectoryToImage<PT>::SerieDataType 
DicomDirectoryToImage<PT>::GetAllDicomSeriesFromSubDir(
            const std::string subdir
        ) 
{
    if (m_Verbose)
    {
        std::cout << " reading series ... " << std::endl;
    }

    SerieDataType allseries; 

    // reading series 
    // uses the first file to obtain series information
    // seriesdescription and seriesid
    typedef GDCMSeriesFileNames                  NamesGeneratorType;
    typename NamesGeneratorType::Pointer namegenerator = NamesGeneratorType::New();
    namegenerator->SetUseSeriesDetails(true);
    namegenerator->AddSeriesRestriction("0008|0021" );
    namegenerator->SetDirectory(subdir);

    const std::vector<std::string> & seriesuid 
        = namegenerator->GetSeriesUIDs();

    // traverse through the series
    for (unsigned int i=0; i<seriesuid.size(); ++i) 
    {

        // read this series
        const std::vector<std::string> filenames 
            = namegenerator->GetFileNames(seriesuid[i].c_str());


        // read only first dicom file of this serie 
        gdcm::ImageReader r;

        r.SetFileName(filenames[0].c_str());
        r.Read();

        const gdcm::DataSet h = r.GetFile().GetDataSet();

        // retrieve seriesdes
        std::string sd("");
        gdcm::Attribute<0x0008,0x103e> atsd;
        atsd.SetFromDataElement(h.GetDataElement(atsd.GetTag()));
        if (!h.GetDataElement(atsd.GetTag()).IsEmpty())
        {
            sd = atsd.GetValue();
        }
        else
        {
            std::cout << "error: reading seriesdescription from file " << filenames[0] << std::endl;
        }
 
        if (m_Verbose)
        {
            std::cout << " [" << i << "]";
            std::cout.width(4);
            std::cout << filenames.size() << " ";
            std::cout.width(40);
            std::cout << sd << "\t" << seriesuid[i].c_str()<< std::endl;
        }
        SD serie;
        serie.nslices = filenames.size();
        serie.description = sd;
        serie.id = seriesuid[i];
        allseries[i] = serie;
    }
    return allseries;
}
// applycriteria
template<class PT>
std::vector<int> DicomDirectoryToImage<PT>
::ApplyCriteria(
            const SerieDataType & sd
        )
{
    if (m_Verbose)
    {
        std::cout << " applying criteria ... " << std::endl; 
        std::cout << " minimum slices      \t" << m_MinimumSlices << std::endl;
        PrintMatchCriterium();
    }
    // given all series, apply the following criteria
    // - match on series description (AND/OR)
    // - apply minimum number of slices

    std::vector<int> sel;

    typename SerieDataType::const_iterator it;

    for (it = sd.begin(); it!=sd.end(); ++it)
    {
        const unsigned int ns = (it->second).nslices;

        std::string de = (it->second).description;
        std::transform(de.begin(), de.end(), de.begin(), tolower);

        if ( ns >= m_MinimumSlices && ns > 1)
        {
            if (m_CriteriumSingle.empty() &&
                    m_CriteriumOr.empty() &&
                    m_CriteriumAnd.empty())
            {
                sel.push_back(it->first);
            }
            else
            {
                // check single criterium
                if (!m_CriteriumSingle.empty())
                {
                    if (de.find(m_CriteriumSingle) != std::string::npos)
                    {
                        sel.push_back(it->first);
                    }
                }
                // or criteria
                bool matchor = false;
                if (!m_CriteriumOr.empty()) 
                {
                    unsigned int i(0);
                    while (!matchor && i<m_CriteriumOr.size())
                    {
                        if (de.find(m_CriteriumOr[i]) != std::string::npos)
                        {
                            matchor = true;
                            sel.push_back(it->first);
                        }
                        ++i;
                    }
                }

                // and 
                if (!matchor && !m_CriteriumAnd.empty())
                {
                    unsigned int i(0);
                    bool matchand = true;
                    while (matchand && i < m_CriteriumAnd.size())
                    {
                        if (de.find(m_CriteriumAnd[i]) == std::string::npos)
                        {
                            matchand = false;
                        }
                        ++i;
                    }
                    if (matchand)
                    {
                        sel.push_back(it->first);
                    }
                }
            }
        } // end check nslices
    }
    if (m_Verbose)
    {
        if (sel.empty())
        {
            std::cout << " selected series     \tnone" << std::endl;
        }
        else
        {
            std::cout << " selected series     \t";
            for (unsigned int i=0 ; i<sel.size(); ++i)
            {
                std::cout << sel[i] << " ";
            }
            std::cout << std::endl;
        }
    }
    return sel;
}
// construct filenames
template<class PT>
std::string DicomDirectoryToImage<PT>
::ConstructBaseName(   const std::string patientname, 
                        const std::string patientid,
                        const std::string acquisitiondate, 
                        const std::string convolutionkernel,
                        const std::string accnumber)
{
    std::string patientnameclean(""); 
    std::string patientidclean("");
    std::string acquisitiondateclean("");
    std::string convolutionkernelclean("");
    std::string accnumberclean("");

    // patientname
    patientnameclean = patientname.substr(0, patientname.find_first_of(","));

    while (patientnameclean.find(" ") != std::string::npos)
        patientnameclean.erase(patientnameclean.find(" "),1);

    while (patientnameclean.find(".") != std::string::npos)
        patientnameclean.erase(patientnameclean.find("."),1);

    while (patientnameclean.find("-") != std::string::npos)
        patientnameclean.erase(patientnameclean.find("-"),1);

    while (patientnameclean.find(",") != std::string::npos)
        patientnameclean.erase(patientnameclean.find(","),1);

    while (patientnameclean.find("^") != std::string::npos)
        patientnameclean.erase(patientnameclean.find("^"),1);

    std::transform (patientnameclean.begin(),
                    patientnameclean.end(),
                    patientnameclean.begin(),
                    tolower);

    // patientid
    patientidclean = patientid;

    while (patientidclean.find(" ") != std::string::npos)
        patientidclean.erase(patientidclean.find(" "),1);

    while (patientidclean.find(".") != std::string::npos)
        patientidclean.erase(patientidclean.find("."),1);

    // acquisitiondate
    acquisitiondateclean = acquisitiondate;

    while (acquisitiondateclean.find(" ") != std::string::npos)
        acquisitiondateclean.erase(acquisitiondateclean.find(" "),1);

    while (acquisitiondateclean.find(".") != std::string::npos)
        acquisitiondateclean.erase(acquisitiondateclean.find("."),1);

    while (acquisitiondateclean.find("-") != std::string::npos)
        acquisitiondateclean.erase(acquisitiondateclean.find("-"),1);

    // convolutionkernel
    convolutionkernelclean = convolutionkernel;

    while (convolutionkernelclean.find(" ") != std::string::npos)
        convolutionkernelclean.erase(convolutionkernelclean.find(" "),1);

    /*
    while (convolutionkernelclean.find("f") != std::string::npos)
        convolutionkernelclean.erase(convolutionkernelclean.find("f"),1);

    while (convolutionkernelclean.find("F") != std::string::npos)
        convolutionkernelclean.erase(convolutionkernelclean.find("F"),1);
        */

    std::transform (convolutionkernelclean.begin(),
                    convolutionkernelclean.end(),
                    convolutionkernelclean.begin(),
                    tolower);

    if (m_UseAccessionNumberInNaming)
    {
        // accnumber
        accnumberclean = accnumber.substr(0, accnumber.find_first_of("."));

        while (accnumberclean.find(" ") != std::string::npos)
            accnumberclean.erase(accnumberclean.find(" "),1);

        while (accnumberclean.find("-") != std::string::npos)
            accnumberclean.erase(accnumberclean.find("-"),1);

        while (accnumberclean.find(",") != std::string::npos)
            accnumberclean.erase(accnumberclean.find(","),1);

        std::transform (accnumberclean.begin(),
                accnumberclean.end(),
                accnumberclean.begin(),
                tolower);
    }

    std::string base(""); 

    // first element
    if (m_Anonymize)
    {
        base = patientidclean;
    }
    else
    {
        base = patientnameclean;
    }

    // patientid
    patientidclean = "_" + patientidclean;
    if (patientidclean != "_" && !m_Anonymize)
    {
        base += patientidclean;
    }

    if (m_UseAccessionNumberInNaming)
    {
        accnumberclean = "_" + accnumberclean;
        if (accnumberclean != "_")
        {
            base += accnumberclean;
        }
    }

    acquisitiondateclean = "_" + acquisitiondateclean;
    if (acquisitiondateclean != "_")
    {
        base += acquisitiondateclean;
    }

    convolutionkernelclean = "_" + convolutionkernelclean;
    if (convolutionkernelclean != "_")
    {
        base += convolutionkernelclean;
    }
    if (base.empty())
    {
        base = "noheaderinfo";
    }

    return base;

}
// writeserie
template<class PT>
bool DicomDirectoryToImage<PT>
::WriteSerie(
            const std::string dir,
            SerieDataType  sd,
            const unsigned int s
        )
{

    std::string msg("");

    std::vector<std::string> emsg;
    bool success(true);

    std::stringstream s1;
    s1 << s; 
    msg = " [" + s1.str() + "]";

    // flowline
    // 1. first retrieving all files for this serie
    //    and do some check on consistency, eg number of files
    // 2. reading first file of serie to retrieve some dcm
    //    tags for constructing file name. This also is used
    //    in 4. to retrieve additional dicom tags
    // 3. importing serie to a volume, writing to volume
    //    using user given extension
    // 4. replacing dcm file (if exist) with one from 2,
    //    modifyid to fit new volume

    // 1. retrieve all files for this series
    typedef GDCMSeriesFileNames   NamesGeneratorType;

    typename NamesGeneratorType::Pointer namegenerator = NamesGeneratorType::New();
    namegenerator->SetUseSeriesDetails(true);
    //namegenerator->AddSeriesRestriction("0008|0021"); // seriesdate
    //namegenerator->AddSeriesRestriction("0028|0010"); // rows
    namegenerator->SetDirectory(dir);

    std::vector<std::string> filenames = 
        namegenerator->GetFileNames(((sd[s]).id).c_str());

    // for sanity do check on dimension for each dcm file in
    // this serie, somehow sometimes dcm files are found with
    // a different x/y dimension. We take the first dcm file
    // as the one with right x/y dimension

    // erases 6th element
    // filenames.erase(filenames.begin()+5);

    if (filenames.size() > 1) // should always be the case
                              // since these series have been selected already,
                              // for sanity we still do a check here
    {
        bool    xyf(false);
        int     X(-1), Y(-1);
        unsigned int    i(0);
        std::vector<unsigned int> erase;
        while (!xyf && i<filenames.size())
        {

            gdcm::ImageReader r;
            r.SetFileName(filenames[i].c_str());
            r.Read();
            const gdcm::DataSet h = r.GetFile().GetDataSet();

            bool a(false), b(false);

            gdcm::Attribute<0x0028,0x0010> atX;
            if (!h.GetDataElement(atX.GetTag()).IsEmpty())
            {
                atX.SetFromDataElement(h.GetDataElement(atX.GetTag()));
                X = atX.GetValue();
                a = true;
            }
            else
            {
                a = false;
            }

            gdcm::Attribute<0x0028,0x0011> atY;
            if (!h.GetDataElement(atY.GetTag()).IsEmpty())
            {
                atY.SetFromDataElement(h.GetDataElement(atY.GetTag()));
                Y = atY.GetValue();
                b = true;
            }
            else
            {
                b = false;
            }

            if (a && b)
            {
                xyf = true;
            }
            else
            {
                // remove this one
                erase.push_back(i);
                // check next
                xyf = false;
                ++i;
            }
        } // end selection first x,y
       
        // check remaining, remove if required
        for (;i<filenames.size(); ++i)
        {
            gdcm::ImageReader r;
            r.SetFileName(filenames[i].c_str());
            r.Read();
            const gdcm::DataSet h = r.GetFile().GetDataSet();

            int x(-1), y(-1);

            gdcm::Attribute<0x0028,0x0010> atX;
            if (!h.GetDataElement(atX.GetTag()).IsEmpty())
            {
                atX.SetFromDataElement(h.GetDataElement(atX.GetTag()));
                x = atX.GetValue();
            }

            gdcm::Attribute<0x0028,0x0011> atY;
            if (!h.GetDataElement(atY.GetTag()).IsEmpty())
            {
                atY.SetFromDataElement(h.GetDataElement(atY.GetTag()));
                y = atY.GetValue();
            }

            if (x!=X || y!=Y)
            {
                erase.push_back(i);
            }
        }
        // do removal
        for (unsigned int i=0; i< erase.size(); ++i)
        {
            filenames.erase(filenames.begin()+(erase[i]));
            emsg.push_back(" WARNING: removal dcmfile because of diff xy size: " + filenames[erase[i]]);  
            success = false;
        }
    }

    std::stringstream s2;
    s2 << filenames.size();
    msg += " " + s2.str() + " ";
 
    // 2. retrieve first gdcm file, this one will be used
    // to retrieve detailed header information, independent of the extension
    // that is selected. Retrieve patientnm, id, acqd, kern
    // intercept and slope
    gdcm::Reader gr;
    gr.SetFileName(filenames[0].c_str());
    gr.Read();

    gdcm::DataSet gh = gr.GetFile().GetDataSet();

    // patientname
    std::string pn("");
    gdcm::Attribute<0x0010,0x0010> atpn;
    if (!gh.GetDataElement(atpn.GetTag()).IsEmpty())
    {
        atpn.SetFromDataElement(gh.GetDataElement(atpn.GetTag()));
        pn = atpn.GetValue();
    }
    else
    {
        emsg.push_back(" ERROR: retrieving patientname from " + filenames[0]);  
        success = false;
    }
    // patientid
    std::string pid("");
    gdcm::Attribute<0x0010,0x0020> atpid;
    if (!gh.GetDataElement(atpid.GetTag()).IsEmpty())
    {
        atpid.SetFromDataElement(gh.GetDataElement(atpid.GetTag()));
        pid = atpid.GetValue();
    }
    else
    {
        emsg.push_back(" ERROR: retrieving patientid from " + filenames[0]);
        success = false;
    }
    // acqdata
    std::string aqd("");
    gdcm::Attribute<0x0008,0x0022> ataqd;
    if (!gh.GetDataElement(ataqd.GetTag()).IsEmpty())
    {
        ataqd.SetFromDataElement(gh.GetDataElement(ataqd.GetTag()));
        aqd = ataqd.GetValue();
    }
    else
    {
        emsg.push_back(" ERROR: retrieving acquistiondate from " + filenames[0]);
        success = false;
    }
    // kernel
    std::string kern("");
    gdcm::Attribute<0x0018,0x1210> atkern;
    if (!gh.GetDataElement(atkern.GetTag()).IsEmpty())
    {
        atkern.SetFromDataElement(gh.GetDataElement(atkern.GetTag()));
        kern = atkern.GetValue();
    }
    else
    {
        emsg.push_back(" ERROR: retrieving kernel from " + filenames[0]);
        success = false;
    }
    // intercept 
    double intercept(NumericTraits<double>::Zero);
    gdcm::Attribute<0x0028,0x1052> atri;
    if (!gh.GetDataElement(atri.GetTag()).IsEmpty())
    {
        atri.SetFromDataElement(gh.GetDataElement(atri.GetTag()));
        intercept = atri.GetValue();
    }
    else
    {
        emsg.push_back(" ERROR: retrieving rescale intercept from " + filenames[0]);
        success = false;
    }

    // slope
    double slope(NumericTraits<double>::One);
    gdcm::Attribute<0x0028,0x1053> atrs;
    if (!gh.GetDataElement(atrs.GetTag()).IsEmpty())
    {
        atrs.SetFromDataElement(gh.GetDataElement(atrs.GetTag()));
        slope = atrs.GetValue();
    }
    else
    {
        emsg.push_back(" ERROR: retrieving rescale slope from " + filenames[0]);
        success = false;
    }
    std::string acc("");
    if (m_UseAccessionNumberInNaming)
    {
        gdcm::Attribute<0x0008,0x0050> atacc;
        if (!gh.GetDataElement(atacc.GetTag()).IsEmpty())
        {
            atacc.SetFromDataElement(gh.GetDataElement(atacc.GetTag()));
            acc = atacc.GetValue();
        }
        else
        {
            emsg.push_back(" ERROR: retrieving accessionnumber from " + filenames[0]);
            success = false;
        }
    }

   // check for presence of dcmname (because this one is always written)
    // if clean is found, then rename the found files with post 000, 
    //    and add count for next
    // if present with postfix, then
    //    get all current files with enumerate postfix, retrieve max
    //    and add for next
    std::string base = ConstructBaseName(pn,pid,aqd,kern,acc);
    struct stat fileinfo;

    const std::string check = m_OutputDirectory + "/" + base + ".dcm";
    const bool basepres( (stat(check.c_str(), &fileinfo) == 0)? true:false);

    const std::string checkp = m_OutputDirectory + "/" + base + "_0000.dcm";
    const bool postpres( (stat(checkp.c_str(), &fileinfo) == 0)? true:false);

    std::string post("");
    if (!basepres && !postpres)
    {
        post = "";
    }
    if (basepres && !postpres)
    {
        post = "_0000";
    }
    if (postpres)
    {
        // retrieve current max nr
        bool f(false);
        int c(0);
        while (!f)
        {
            const int n(c+1);
            std::stringstream sc,sn;
            sc << c;
            sn << n;
            std::string current(sc.str());
            std::string next(sn.str());
            if (current.size() == 1)
            {
                current = "000" + current;
            }
            if (current.size() == 2)
            {
                current = "00" + current;
            }
            if (current.size() == 3)
            {
                current = "0" + current;
            }
            current = "_" + current;
            if (next.size() == 1)
            {
                next = "000" + next;
            }
            if (next.size() == 2)
            {
                next = "00" + next;
            }
            if (next.size() == 3)
            {
                next = "0" + next;
            }
            next = "_" + next;

            const std::string filec(m_OutputDirectory + "/" + base + current + ".dcm");
            const std::string filen(m_OutputDirectory + "/" + base + next + ".dcm");

            const bool cp( (stat(filec.c_str(), &fileinfo) == 0)? true:false);
            const bool np( (stat(filen.c_str(), &fileinfo) == 0)? true:false);

            if ( (cp && !np) || c == 9998)
            {
                f = true;
            }
            else
            {
                c++;
            }
        }

        c++;
        std::stringstream sp;
        sp << c;
        post = sp.str();
        if (post.size() == 1)
        {
            post = "000" + post;
        }
        if (post.size() == 2)
        {
            post = "00" + post;
        }
        if (post.size() == 3)
        {
            post = "0" + post;
        }
        post = "_" + post;
    }

    // construct names
    base += post;
    const std::string volname = base + "." + m_FileNameExtension;
    const std::string dcmname = base + ".dcm";

    // 3. VOLUME
    //
    // actually read series/and writing
    // If the extension is selected as mevis dcm/tiff
    // then it automatically writes a dcm file, which needs
    // to be replaced later on
    typedef GDCMImageIO   ImageIOType;
    typename ImageIOType::Pointer  dicomio = ImageIOType::New();

    // somehow dicom rescaling is always applied!
    // and reflected in the corresponding dcm file
    // therefore we always use SHORT pixeltype, and
    // adjust rescale inter/slope in dcm file

    typedef Image<PT,3>             IT;
    typedef ImageSeriesReader<IT>   RT;

    typename RT::Pointer r = RT::New();
    r->SetImageIO(dicomio);
    r->SetFileNames(filenames);
    try {
        r->Update();
    }
    catch (ExceptionObject &excp)
    {
        emsg.push_back(" ERROR: reading filesserie, exception caught ");
        std::stringstream s;
        s << excp;
        emsg.push_back(" " + s.str());
        success = false;
    }

    std::stringstream s3,s4;
    s3 << intercept;
    s4 << slope;
    msg += "interc/slope " + s3.str() + " " + s4.str() + " auto applied, ";
    msg += "writing: " + volname + ", ";

    typedef ImageFileWriter<IT> WT;
    typename WT::Pointer w = WT::New();
    w->SetFileName(m_OutputDirectory + "/" + volname);
    w->SetInput(r->GetOutput());
    try {
        w->Update();
    }
    catch (ExceptionObject &excp)
    {
        emsg.push_back(" ERROR: writing volume, exception caught ");
        std::stringstream s;
        s << excp;
        emsg.push_back(" " + s.str());
        success = false;
    }

    // 4. DCM
    //
    // if a dcm file exists (because in 3 extension dcm/tif has been
    // selected) then this one is read and extended with the list
    // of dcm tags below. Note: gdcm::ImageReader below does even work 
    // if a dcm file does not exist, so no explicit check on dcm file
    // existence

    gdcm::ImageReader newdcm; 
    newdcm.SetFileName((m_OutputDirectory + "/" + dcmname).c_str());
    newdcm.Read();
    gdcm::DataSet newdcmdataset = newdcm.GetFile().GetDataSet();

    // copied from original dcm to newdcm
    gdcm::Attribute<0x0008,0x0008> atimagetype;
    if (!gh.GetDataElement(atimagetype.GetTag()).IsEmpty())
    {
        atimagetype.SetFromDataElement(gh.GetDataElement(atimagetype.GetTag()));
        newdcmdataset.Replace(atimagetype.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0008 (image type) from dcm-file:";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0016> atsopclassuid;
    if (!gh.GetDataElement(atsopclassuid.GetTag()).IsEmpty())
    {
        atsopclassuid.SetFromDataElement(gh.GetDataElement(atsopclassuid.GetTag()));
        newdcmdataset.Replace(atsopclassuid.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0016 (sop class uid) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0018> atsopinstanceuid;
    if (!gh.GetDataElement(atsopinstanceuid.GetTag()).IsEmpty())
    {
        atsopinstanceuid.SetFromDataElement(gh.GetDataElement(atsopinstanceuid.GetTag()));
        newdcmdataset.Replace(atsopinstanceuid.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0018 (sop instance uid) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0020> atstudydate;
    if (!gh.GetDataElement(atstudydate.GetTag()).IsEmpty())
    {
        atstudydate.SetFromDataElement(gh.GetDataElement(atstudydate.GetTag()));
        newdcmdataset.Replace(atstudydate.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0020 (study date) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0021> atseriedate;
    if (!gh.GetDataElement(atseriedate.GetTag()).IsEmpty())
    {
        atseriedate.SetFromDataElement(gh.GetDataElement(atseriedate.GetTag()));
        newdcmdataset.Replace(atseriedate.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0021 (serie date) from dcm-file"; 
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0022> atacqdate;
    if (!gh.GetDataElement(atacqdate.GetTag()).IsEmpty())
    {
        atacqdate.SetFromDataElement(gh.GetDataElement(atacqdate.GetTag()));
        newdcmdataset.Replace(atacqdate.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0022 (acq date) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0023> atcontdate;
    if (!gh.GetDataElement(atcontdate.GetTag()).IsEmpty())
    {
        atcontdate.SetFromDataElement(gh.GetDataElement(atcontdate.GetTag()));
        newdcmdataset.Replace(atcontdate.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0022 (contdate) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0008,0x0030> atstudytime;
    if (!gh.GetDataElement(atstudytime.GetTag()).IsEmpty())
    {
        atstudytime.SetFromDataElement(gh.GetDataElement(atstudytime.GetTag()));
        newdcmdataset.Replace(atstudytime.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0030 (studytime) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0008,0x0031> atserietime;
    if (!gh.GetDataElement(atserietime.GetTag()).IsEmpty())
    {
        atserietime.SetFromDataElement(gh.GetDataElement(atserietime.GetTag()));
        newdcmdataset.Replace(atserietime.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0031 (serietime) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0008,0x0032> atacqtime;
    if (!gh.GetDataElement(atacqtime.GetTag()).IsEmpty())
    {
        atacqtime.SetFromDataElement(gh.GetDataElement(atacqtime.GetTag()));
        newdcmdataset.Replace(atacqtime.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0032 (acqtime) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0008,0x0033> atconttime;
    if (!gh.GetDataElement(atconttime.GetTag()).IsEmpty())
    {
        atconttime.SetFromDataElement(gh.GetDataElement(atconttime.GetTag()));
        newdcmdataset.Replace(atconttime.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0033 (conttime) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0050> ataccnr;
    if (!gh.GetDataElement(ataccnr.GetTag()).IsEmpty())
    {
        ataccnr.SetFromDataElement(gh.GetDataElement(ataccnr.GetTag()));
        newdcmdataset.Replace(ataccnr.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0050 (accession nr) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0008,0x0060> atmod;
    if (!gh.GetDataElement(atmod.GetTag()).IsEmpty())
    {
        atmod.SetFromDataElement(gh.GetDataElement(atmod.GetTag()));
        newdcmdataset.Replace(atmod.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0060 (modality) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x0070> atman;
    if (!gh.GetDataElement(atman.GetTag()).IsEmpty())
    {
        atman.SetFromDataElement(gh.GetDataElement(atman.GetTag()));
        newdcmdataset.Replace(atman.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,0070 (manufacturer) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0008,0x0090> atref;
    if (!gh.GetDataElement(atref.GetTag()).IsEmpty())
    {
        atref.SetFromDataElement(gh.GetDataElement(atref.GetTag()));
        newdcmdataset.Replace(atref.GetAsDataElement());
    }    
    else
    {
        std::cout << "dicomimport: error reading 0008,0090 (ref) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x1010> atstat;
    if (!gh.GetDataElement(atstat.GetTag()).IsEmpty())
    {
        atstat.SetFromDataElement(gh.GetDataElement(atstat.GetTag()));
        newdcmdataset.Replace(atstat.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,1010 (stat) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x1030> atstudydescr;
    if (!gh.GetDataElement(atstudydescr.GetTag()).IsEmpty())
    {
        atstudydescr.SetFromDataElement(gh.GetDataElement(atstudydescr.GetTag()));
        newdcmdataset.Replace(atstudydescr.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,1030 (studydescr) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x103e> atseriedescr;
    if (!gh.GetDataElement(atseriedescr.GetTag()).IsEmpty())
    {
        atseriedescr.SetFromDataElement(gh.GetDataElement(atseriedescr.GetTag()));
        newdcmdataset.Replace(atseriedescr.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,103e (seriedescr) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0008,0x1048> atphys;
    if (!gh.GetDataElement(atphys.GetTag()).IsEmpty())
    {
        atphys.SetFromDataElement(gh.GetDataElement(atphys.GetTag()));
        newdcmdataset.Replace(atphys.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,1048 (phys) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0008,0x1090> atmodel;
    if (!gh.GetDataElement(atmodel.GetTag()).IsEmpty())
    {
        atmodel.SetFromDataElement(gh.GetDataElement(atmodel.GetTag()));
        newdcmdataset.Replace(atmodel.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0008,1090 (model) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0010,0x0020> atpatientid;
    if (!gh.GetDataElement(atpatientid.GetTag()).IsEmpty())
    {
        atpatientid.SetFromDataElement(gh.GetDataElement(atpatientid.GetTag()));
        newdcmdataset.Replace(atpatientid.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0010,0020 (patientid) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x0015> atbody;
    if (!gh.GetDataElement(atbody.GetTag()).IsEmpty())
    {
        atbody.SetFromDataElement(gh.GetDataElement(atbody.GetTag()));
        newdcmdataset.Replace(atbody.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,0015 (body) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x0060> atkvp;
    if (!gh.GetDataElement(atkvp.GetTag()).IsEmpty())
    {
        atkvp.SetFromDataElement(gh.GetDataElement(atkvp.GetTag()));
        newdcmdataset.Replace(atkvp.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,0060 (kvp) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x0090> atcoldiam;
    if (!gh.GetDataElement(atcoldiam.GetTag()).IsEmpty())
    {
        atcoldiam.SetFromDataElement(gh.GetDataElement(atcoldiam.GetTag()));
        newdcmdataset.Replace(atcoldiam.GetAsDataElement());
    }    
    else
    {
        std::cout << "dicomimport: error reading 0018,0090 (coldiam) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1000> atdevicesn;
    if (!gh.GetDataElement(atdevicesn.GetTag()).IsEmpty())
    {
        atdevicesn.SetFromDataElement(gh.GetDataElement(atdevicesn.GetTag()));
        newdcmdataset.Replace(atdevicesn.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1000 (devicesn) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1020> atsoftwarevs;
    if (!gh.GetDataElement(atsoftwarevs.GetTag()).IsEmpty())
    {
        atsoftwarevs.SetFromDataElement(gh.GetDataElement(atsoftwarevs.GetTag()));
        newdcmdataset.Replace(atsoftwarevs.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1020 (softwarevs) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1030> atprotocol;
    if (!gh.GetDataElement(atprotocol.GetTag()).IsEmpty())
    {
        atprotocol.SetFromDataElement(gh.GetDataElement(atprotocol.GetTag()));
        newdcmdataset.Replace(atprotocol.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1030 (protocol) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1100> atrecdiam;
    if (!gh.GetDataElement(atrecdiam.GetTag()).IsEmpty())
    {
        atrecdiam.SetFromDataElement(gh.GetDataElement(atrecdiam.GetTag()));
        newdcmdataset.Replace(atrecdiam.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1100 (recdiam) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1110> atsourcedet;
    if (!gh.GetDataElement(atsourcedet.GetTag()).IsEmpty())
    {
        atsourcedet.SetFromDataElement(gh.GetDataElement(atsourcedet.GetTag()));
        newdcmdataset.Replace(atsourcedet.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1110 (sourcedet) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1111> atsourcepat;
    if (!gh.GetDataElement(atsourcepat.GetTag()).IsEmpty())
    {
        atsourcepat.SetFromDataElement(gh.GetDataElement(atsourcepat.GetTag()));
        newdcmdataset.Replace(atsourcepat.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1111 (sourcepat) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1120> atgantrytilt;
    if (!gh.GetDataElement(atgantrytilt.GetTag()).IsEmpty())
    {
        atgantrytilt.SetFromDataElement(gh.GetDataElement(atgantrytilt.GetTag()));
        newdcmdataset.Replace(atgantrytilt.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1120 (sourcepat) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1130> attableh;
    if (!gh.GetDataElement(attableh.GetTag()).IsEmpty())
    {
        attableh.SetFromDataElement(gh.GetDataElement(attableh.GetTag()));
        newdcmdataset.Replace(attableh.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1130 (tableheight) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1140> atrotdir;
    if (!gh.GetDataElement(atrotdir.GetTag()).IsEmpty())
    {
        atrotdir.SetFromDataElement(gh.GetDataElement(atrotdir.GetTag()));
        newdcmdataset.Replace(atrotdir.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1140 (rot dir) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0018,0x1150> atexptime;
    if (!gh.GetDataElement(atexptime.GetTag()).IsEmpty())
    {
        atexptime.SetFromDataElement(gh.GetDataElement(atexptime.GetTag()));
        newdcmdataset.Replace(atexptime.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1150 (exp time) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1151> atmas;
    if (!gh.GetDataElement(atmas.GetTag()).IsEmpty())
    {
        atmas.SetFromDataElement(gh.GetDataElement(atmas.GetTag()));
        newdcmdataset.Replace(atmas.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1151 (mas) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0018,0x1152> atexp;
    if (!gh.GetDataElement(atexp.GetTag()).IsEmpty())
    {
        atexp.SetFromDataElement(gh.GetDataElement(atexp.GetTag()));
        newdcmdataset.Replace(atexp.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1152 (exp) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0018,0x1160> atfilter;
    if (!gh.GetDataElement(atfilter.GetTag()).IsEmpty())
    {
        atfilter.SetFromDataElement(gh.GetDataElement(atfilter.GetTag()));
        newdcmdataset.Replace(atfilter.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1160 (filter) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1170> atgenpow;
    if (!gh.GetDataElement(atgenpow.GetTag()).IsEmpty())
    {
        atgenpow.SetFromDataElement(gh.GetDataElement(atgenpow.GetTag()));
        newdcmdataset.Replace(atgenpow.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1170 (genpow) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1190> atfocal;
    if (!gh.GetDataElement(atfocal.GetTag()).IsEmpty())
    {
        atfocal.SetFromDataElement(gh.GetDataElement(atfocal.GetTag()));
        newdcmdataset.Replace(atfocal.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1190 (focal) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0018,0x1200> atcaldate;
    if (!gh.GetDataElement(atcaldate.GetTag()).IsEmpty())
    {
        atcaldate.SetFromDataElement(gh.GetDataElement(atcaldate.GetTag()));
        newdcmdataset.Replace(atcaldate.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1200 (caldate) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0018,0x1201> atcaltime;
    if (!gh.GetDataElement(atcaltime.GetTag()).IsEmpty())
    {
        atcaltime.SetFromDataElement(gh.GetDataElement(atcaltime.GetTag()));
        newdcmdataset.Replace(atcaltime.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1201 (caltime) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0018,0x1210> atconvkern;
    if (!gh.GetDataElement(atconvkern.GetTag()).IsEmpty())
    {
        atconvkern.SetFromDataElement(gh.GetDataElement(atconvkern.GetTag()));
        newdcmdataset.Replace(atconvkern.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,1210 (convkern) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0018,0x5100> atpatpos;
    if (!gh.GetDataElement(atpatpos.GetTag()).IsEmpty())
    {
        atpatpos.SetFromDataElement(gh.GetDataElement(atpatpos.GetTag()));
        newdcmdataset.Replace(atpatpos.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0018,5100 (patpos) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0020,0x000d> atstudyuid;
    if (!gh.GetDataElement(atstudyuid.GetTag()).IsEmpty())
    {
        atstudyuid.SetFromDataElement(gh.GetDataElement(atstudyuid.GetTag()));
        newdcmdataset.Replace(atstudyuid.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0020,000d (studyuid) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0020,0x000e> atserieuid;
    if (!gh.GetDataElement(atserieuid.GetTag()).IsEmpty())
    {
        atserieuid.SetFromDataElement(gh.GetDataElement(atserieuid.GetTag()));
        newdcmdataset.Replace(atserieuid.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0020,000e (serieuid) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0020,0x0010> atstudyid;
    if (!gh.GetDataElement(atstudyid.GetTag()).IsEmpty())
    {
        atstudyid.SetFromDataElement(gh.GetDataElement(atstudyid.GetTag()));
        newdcmdataset.Replace(atstudyid.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0020,0010 (studyid) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0020,0x0011> atseriesnr;
    if (!gh.GetDataElement(atseriesnr.GetTag()).IsEmpty())
    {
        atseriesnr.SetFromDataElement(gh.GetDataElement(atseriesnr.GetTag()));
        newdcmdataset.Replace(atseriesnr.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0020,0011 (seriesnr) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0020,0x0012> atacqnr;
    if (!gh.GetDataElement(atacqnr.GetTag()).IsEmpty())
    {
        atacqnr.SetFromDataElement(gh.GetDataElement(atacqnr.GetTag()));
        newdcmdataset.Replace(atacqnr.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0020,0012 (acqnr) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0020,0x0105> attemp;
    if (!gh.GetDataElement(attemp.GetTag()).IsEmpty())
    {
        attemp.SetFromDataElement(gh.GetDataElement(attemp.GetTag()));
        newdcmdataset.Replace(attemp.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0020,0105 (temp) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0028,0x1050> atwindc;
    if (!gh.GetDataElement(atwindc.GetTag()).IsEmpty())
    {
        atwindc.SetFromDataElement(gh.GetDataElement(atwindc.GetTag()));
        newdcmdataset.Replace(atwindc.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0028,1050 (windc) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    gdcm::Attribute<0x0028,0x1051> atwindw;
    if (!gh.GetDataElement(atwindw.GetTag()).IsEmpty())
    {
        atwindw.SetFromDataElement(gh.GetDataElement(atwindw.GetTag()));
        newdcmdataset.Replace(atwindw.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0028,1051 (windw) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }

    gdcm::Attribute<0x0028,0x1055> atwind;
    if (!gh.GetDataElement(atwind.GetTag()).IsEmpty())
    {
        atwind.SetFromDataElement(gh.GetDataElement(atwind.GetTag()));
        newdcmdataset.Replace(atwind.GetAsDataElement());
    }
    else
    {
        std::cout << "dicomimport: error reading 0028,1055 (wind) from dcm-file";
        std::cout << "\t" + filenames[0] << std::endl;
    }


    if (!m_Anonymize)
    {
        gdcm::Attribute<0x0008,0x0080> atinst;
        if (!gh.GetDataElement(atinst.GetTag()).IsEmpty())
        {
            atinst.SetFromDataElement(gh.GetDataElement(atinst.GetTag()));
            newdcmdataset.Replace(atinst.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0008,0080 (inst) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }


        gdcm::Attribute<0x0008,0x0081> ataddr;
        if (!gh.GetDataElement(ataddr.GetTag()).IsEmpty())
        {
            ataddr.SetFromDataElement(gh.GetDataElement(ataddr.GetTag()));
            newdcmdataset.Replace(ataddr.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0008,0081 (addr) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }

        gdcm::Attribute<0x0010,0x0010> atpn;
        if (!gh.GetDataElement(atpn.GetTag()).IsEmpty())
        {
            atpn.SetFromDataElement(gh.GetDataElement(atpn.GetTag()));
            newdcmdataset.Replace(atpn.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0010,0010 (pn) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }

        gdcm::Attribute<0x0010,0x0030> atpb;
        if (!gh.GetDataElement(atpb.GetTag()).IsEmpty())
        {
            atpb.SetFromDataElement(gh.GetDataElement(atpb.GetTag()));
            newdcmdataset.Replace(atpb.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0010,0030 (pb) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }

        gdcm::Attribute<0x0010,0x0040> atps;
        if (!gh.GetDataElement(atps.GetTag()).IsEmpty())
        {
            atps.SetFromDataElement(gh.GetDataElement(atps.GetTag()));
            newdcmdataset.Replace(atps.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0010,0040 (ps) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }

        gdcm::Attribute<0x0010,0x1010> atpa;
        if (!gh.GetDataElement(atpa.GetTag()).IsEmpty())
        {
            atpa.SetFromDataElement(gh.GetDataElement(atpa.GetTag()));
            newdcmdataset.Replace(atpa.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0010,1010 (pa) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }

        gdcm::Attribute<0x0010,0x1040> atpad;
        if (!gh.GetDataElement(atpad.GetTag()).IsEmpty())
        {
            atpad.SetFromDataElement(gh.GetDataElement(atpad.GetTag()));
            newdcmdataset.Replace(atpad.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0010,1040 (pad) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }

        gdcm::Attribute<0x0010,0x1060> atmother;
        if (!gh.GetDataElement(atmother.GetTag()).IsEmpty())
        {
            atmother.SetFromDataElement(gh.GetDataElement(atmother.GetTag()));
            newdcmdataset.Replace(atmother.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0010,1060 (patientsmothersbirthname) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }

        gdcm::Attribute<0x0010,0x21c0> atpp;
        if (!gh.GetDataElement(atpp.GetTag()).IsEmpty())
        {
            atpp.SetFromDataElement(gh.GetDataElement(atpp.GetTag()));
            newdcmdataset.Replace(atpp.GetAsDataElement());
        }
        else
        {
            std::cout << "dicomimport: error reading 0010,21c0(pp) from dcm-file";
            std::cout << "\t" + filenames[0] << std::endl;
        }
    }
 
    // always replaced
    // comments
    gdcm::Attribute<0x0020,0x4000> atc;
    const std::string v(Version::GetITKVersion());
    const std::string c = "ITK " + v + " DicomDirectoryToImage";
    atc.SetValue(c);
    newdcmdataset.Replace(atc.GetAsDataElement());

    // intercept/slope
    // Since dicom rescaling is automatically 
    // applied we always adjust values in dcm
    // file
    gdcm::Attribute<0x0028,0x1052> atwri;
    gdcm::Attribute<0x0028,0x1053> atwrs;
    atwri.SetValue(NumericTraits<double>::Zero);
    atwrs.SetValue(NumericTraits<double>::One);

    newdcmdataset.Replace(atwri.GetAsDataElement());
    newdcmdataset.Replace(atwrs.GetAsDataElement());

    // write (overwites dcm file if already exist, because
    // dcm/tiff extension has been selected in step 3.)
    msg += dcmname + "\n";
    gdcm::Writer gw;
    gw.SetCheckFileMetaInformation(false);
    gw.GetFile().SetDataSet(newdcmdataset);
    gw.SetFileName((m_OutputDirectory + "/" + dcmname).c_str());
    gw.Write();

    if (!success)
    {
        for (unsigned int i=0; i<emsg.size(); ++i)
        {
            std::cout << emsg[i] << std::endl;
        }
    }
    m_LogMessages.push_back(msg);
    if (m_Verbose)
    {
        std::cout << msg;
    }

    return success;
}
// getsubdirs
template<class PT>
void DicomDirectoryToImage<PT>
::GetSubDirs(
        const std::string        inputdir)
{    
    // check ifdir
    if (!itksys::SystemTools::FileIsDirectory(
                        itksys::SystemTools::GetRealPath(inputdir.c_str()).c_str()))
    {
        return;
    }
 
    Directory::Pointer dir = Directory::New(); 
    dir->Load(inputdir.c_str());

    unsigned int nf(0);

    for (unsigned int i=0; i<dir->GetNumberOfFiles(); ++i)
    {
        const std::string cf(dir->GetFile(i));

        if (cf != "." && cf != "..")
        {
            const std::string d(inputdir + "/" + cf);

            if (itksys::SystemTools::FileIsDirectory(
                        itksys::SystemTools::GetRealPath(d.c_str()).c_str()))
            {
                // is sub dir
                if (m_Recursive)
                {
                    GetSubDirs(d);
                }
            }
            else
            {
                // is file
                nf++;
            }
        }
    }
    if (nf > 1)
    {
        m_DirsToImport.push_back(inputdir);
    }
}
 
//main: processdirectory
template<class PT>
bool DicomDirectoryToImage<PT>
::ProcessDirectory()
{    
    if (m_Verbose) 
    {
        std::cout << std::endl << "begin dicomdirectorytoimage ..." << std::endl;
    }
    TimeProbe clock;
    clock.Start();
 
    // some cleaning and testing of input vars
    if (m_FileNameExtension.empty())
    {
        m_FileNameExtension = "tif";
    }
    std::transform (m_FileNameExtension.begin(),
            m_FileNameExtension.end(),
            m_FileNameExtension.begin(),
            tolower);

    while (m_FileNameExtension.find(".") != std::string::npos)
    {
        m_FileNameExtension.erase(m_FileNameExtension.find("."),1);
    }

    if (m_DicomDirectory.empty())
    {
        std::cout << "error: no inputdirectory given" << std::endl;
        if (m_Verbose)
        {
            std::cout << "end dicomdirectorytoimage" << std::endl;
        }
        return false;
    }
    if (m_OutputDirectory.empty())
    {
        std::cout << "error: no output directory given" << std::endl;
        if (m_Verbose)
        {
            std::cout << "end dicomdirectorytoimage" << std::endl;
        }
        return false;
    }
    if (!itksys::SystemTools::MakeDirectory(m_OutputDirectory.c_str()))
    {
        std::cout << "error: creating output directory" << std::endl;
        if (m_Verbose)
        {
            std::cout << "end dicomdirectorytoimage" << std::endl;
        }
        return false;
    }

    if (m_Verbose)
    {
        std::cout << "inputdirectory      \t" << m_DicomDirectory << std::endl;
        std::cout << "recursive           \t" << m_Recursive << std::endl;
        std::cout << "outputdirectory     \t" << m_OutputDirectory << std::endl;

        if (!m_MatchOnSeries.empty())
        {
            std::cout << "criteria series     \t" << m_MatchOnSeries<< std::endl;
        }
        std::cout << "minimal nslices     \t" << m_MinimumSlices << std::endl;
        std::cout << "filename extension  \t" << m_FileNameExtension<< std::endl;
        std::cout << "accession in naming \t" << m_UseAccessionNumberInNaming<< std::endl;
    }

 
    // interpret match criterium
   InterpretMatchCriterium();

    // Extract all subdirectories: only add those subdirectories to be
    // imported if it contains at least two files
    // First check ifdir
    if ( !itksys::SystemTools::FileIsDirectory( 
                itksys::SystemTools::GetRealPath(m_DicomDirectory.c_str()).c_str())
            )
    {
        std::cout << "error: " << m_DicomDirectory << " is not a directory" << std::endl;
        if (m_Verbose)
        {
            std::cout << "end dicomdirectorytoimage" << std::endl;
        }
        return false;
    }

    // recursively traverse all directories (note
    // symbolic links are included), and only include
    // directories containing at least two files
    m_DirsToImport.clear();
    GetSubDirs(m_DicomDirectory);

    if (m_Verbose)
    {
        std::cout << "totaldirs to import \t" << m_DirsToImport.size() << std::endl;
    }

    // iterate over dirs
    if (m_Verbose)
    {
        std::cout << "start ... " << std::endl;
    }
 
    std::sort(m_DirsToImport.begin(), m_DirsToImport.end());

    for (unsigned int i=0; i<m_DirsToImport.size(); ++i)
    {
        const std::string cd(m_DirsToImport[i]);

        if (m_Verbose)
        {
            std::cout << std::endl;
            std::cout << static_cast<int>(100.0 * static_cast<float>(i+1))
                            /static_cast<float>(m_DirsToImport.size()) << "% ";
            std::cout << " processing directory " << cd << " ... " << std::endl;
        }

        // retreive all series
        const SerieDataType sd = GetAllDicomSeriesFromSubDir(cd);
        // select series based on criteria
        const std::vector<int> sel = ApplyCriteria(sd);

        // writing series
        if (!sel.empty())
        {
            if (m_Verbose)
            {
                std::cout << " writing series ... " << std::endl;
            }
            for (unsigned int i=0; i< sel.size(); ++i)
            {
                WriteSerie(cd,sd,sel[i]);
            }
        }
    }

    std::cout << "Summary: " << std::endl; 
    for (unsigned int i=0; i<m_LogMessages.size(); ++i)
    {
        std::cout << m_LogMessages[i];
    }
    std::cout << "total data converted\t" << m_LogMessages.size() << std::endl;

    clock.Stop();

    if (m_Verbose)
    {
        std::cout << "end dicomdirectorytoimage" << std::endl;
        std::cout << "elapsed time        \t" << clock.GetMeanTime() << std::endl;
    }
    return true;
}
 

} // end namespace itk

#endif
