/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: dicomimport.cxx,v $
  Language:  C++
  Date:      $Date: 2009/11/04 09:47:41 $
  Version:   $Revision: 1.3 $

=========================================================================*/
#if defined (_MSC_VER)
#pragma warning (disable: 4786)
#endif
#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


#include "itkDicomDirectoryToImage.h"
#include <iostream>


int main (int argc, char ** argv)
{
    itk::DicomDirectoryToImage<short int> import(argv[1]);
    import.SetVerbose(true);
    import.SetRecursive(true);
    import.SetOutputDirectory(argv[2]);
    import.SetMatchOnSeries("CTA");
    import.SetMinimumSlices(10);
    import.SetAnonymize(false);
    import.SetUseAccessionNumberInNaming(true);
    import.SetFileNameExtension("mhd");
    const bool s = import.ProcessDirectory();

    return s;
}

