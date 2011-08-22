/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDicomDirectoryToImage.h,v $
  Language:  C++
  Date:      $Date: 2009/11/01 21:06:52 $
  Version:   $Revision: 1.3 $
=========================================================================*/
#ifndef __itkDicomDirectoryToImage_h
#define __itkDicomDirectoryToImage_h

#include <string>
#include <vector>
#include <map>

namespace itk
{

/** \class DicomDirectoryToImage
 *
 * This class takes an input directory consisting of (raw) dicom files
 * and converts all series to volumetric data (3D only for the moment)
 * and applies consisting filenaming based on the dicom header information.
 *
 * - recursively traversing inputdir, only including those directories
 *   containing at least two files
 * - for each conversion we alway write one dcm file (which
 *   is taken from the first file in the serie)
 * - no files are overwritten, but postfixed with number
 * - templated over output pixel type
 * - requires gdcm 2.x
 * - intercept/slope automatically applied
 * - filenameing,using pidtoctatia number mapping
 * - before distribution with IJ, remove code belonging to ctatia
 *   number mapping
 *
 * nov 1,2009
 *    Release 1.0 for IJ (removed boost dependency)
 *
 *
 */


template<class PT = short int>
class DicomDirectoryToImage
{

public:

    DicomDirectoryToImage(const std::string dirname);
    ~DicomDirectoryToImage();

    // main
    bool ProcessDirectory();

    // helpers
    void SetOutputDirectory(const std::string &a)
    {
        m_OutputDirectory = a;
    }
    void SetVerbose (const bool & a)
    { 
        m_Verbose = a; 
    }
    void SetMatchOnSeries(const std::string &a)
    {
        m_MatchOnSeries = a;
    }
    void SetMinimumSlices(const unsigned int &a)
    {
        m_MinimumSlices = a;
    }
    void SetAnonymize (const bool &a)
    {
        m_Anonymize = a;
    }
    void SetCSVMappingFile(const std::string &a)
    {
        m_CSVMappingFile = a;
    }
    void SetUseAccessionNumberInNaming(const bool & a)
    { 
        m_UseAccessionNumberInNaming = a; 
    }
    void SetFileNameExtension(const std::string & a)
    { 
        m_FileNameExtension = a; 
    }
    void SetRecursive(const bool & a)
    { 
        m_Recursive = a; 
    }
 
 
private:

  bool                          m_Verbose;
  std::string                   m_DicomDirectory;
  std::string                   m_OutputDirectory;
  std::string                   m_MatchOnSeries;
  std::string                   m_InterpretedCriteria;
  unsigned int                  m_MinimumSlices;
  bool                          m_Anonymize;
  std::string                   m_CSVMappingFile;
  bool                          m_UseAccessionNumberInNaming;
  std::string                   m_FileNameExtension;

  std::string                   m_CriteriumSingle;
  std::vector<std::string>      m_CriteriumAnd;
  std::vector<std::string>      m_CriteriumOr;
  std::vector<std::string>      m_DirsToImport;
  bool                          m_Recursive;

  std::vector<std::string>      m_LogMessages;

  // seriesdatatype
  struct SD
  {
      int         nslices;
      std::string description;
      std::string id;
  };
  typedef std::map<int,SD> SerieDataType;

 
  // internal functions
  bool InterpretMatchCriterium();

  void PrintMatchCriterium();

  SerieDataType  GetAllDicomSeriesFromSubDir( const std::string);

  std::vector<int> ApplyCriteria( const SerieDataType &);

  std::string ConstructBaseName( 
          const std::string, 
          const std::string, 
          const std::string, 
          const std::string,
          const std::string);

  bool WriteSerie (
          const std::string, 
          SerieDataType , 
          const unsigned int);

  void GetSubDirs(
          const std::string);
 
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDicomDirectoryToImage.txx"
#endif

#endif
