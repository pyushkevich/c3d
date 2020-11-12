/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertImageND.cxx
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/
#include "ConvertAPI.h"
#include "ConvertImageND.h"
#include "itkSmartPointer.h"
#include "itkImage.h"
#include <cstdarg>

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#else
#include <wordexp.h>
#endif

/** 
 * Code to split a string into argc/argv with some shell expansion
 * from http://stackoverflow.com/questions/1706551/parse-string-into-argv-argc
 */
char **split_commandline(const char *cmdline, int *argc)
{
  int i;
  char **argv = NULL;
  assert(argc);

  if (!cmdline)
    {
    return NULL;
    }

  // Posix.
#ifndef _WIN32
    {
    wordexp_t p;

    // Note! This expands shell variables.
    if (wordexp(cmdline, &p, 0))
      {
      return NULL;
      }

    *argc = p.we_wordc;

    if (!(argv = (char **) calloc(*argc, sizeof(char *))))
      {
      goto fail;
      }

    for (i = 0; i < p.we_wordc; i++)
      {
      if (!(argv[i] = strdup(p.we_wordv[i])))
        {
        goto fail;
        }
      }

    wordfree(&p);

    return argv;
fail:
    wordfree(&p);
    }
#else // WIN32
    {
    wchar_t **wargs = NULL;
    size_t needed = 0;
    wchar_t *cmdlinew = NULL;
    size_t len = strlen(cmdline) + 1;

    if (!(cmdlinew = (wchar_t *) calloc(len, sizeof(wchar_t))))
      goto fail;

    if (!MultiByteToWideChar(CP_ACP, 0, cmdline, -1, cmdlinew, len))
      goto fail;

    if (!(wargs = CommandLineToArgvW(cmdlinew, argc)))
      goto fail;

    if (!(argv = (char **) calloc(*argc, sizeof(char *))))
      goto fail;

    // Convert from wchar_t * to ANSI char *
    for (i = 0; i < *argc; i++)
      {
      // Get the size needed for the target buffer.
      // CP_ACP = Ansi Codepage.
      needed = WideCharToMultiByte(CP_ACP, 0, wargs[i], -1,
        NULL, 0, NULL, NULL);

      if (!(argv[i] = (char *) malloc(needed)))
        goto fail;

      // Do the conversion.
      needed = WideCharToMultiByte(CP_ACP, 0, wargs[i], -1,
        argv[i], needed, NULL, NULL);
      }

    if (wargs) LocalFree(wargs);
    if (cmdlinew) free(cmdlinew);
    return argv;

fail:
    if (wargs) LocalFree(wargs);
    if (cmdlinew) free(cmdlinew);
    }
#endif // WIN32

  if (argv)
    {
    for (i = 0; i < *argc; i++)
      {
      if (argv[i])
        {
        free(argv[i]);
        }
      }

    free(argv);
    }

  return NULL;
}

template <class TPixel, unsigned int VDim>
ConvertAPI<TPixel, VDim>
::ConvertAPI()
{
  m_Converter = new ConverterType();
}

template <class TPixel, unsigned int VDim>
ConvertAPI<TPixel,VDim>
::~ConvertAPI()
{
  delete m_Converter;
}

template <class TPixel, unsigned int VDim>
void 
ConvertAPI<TPixel,VDim>
::AddImage(const char *varname, ImageType *image)
{
  typedef typename ConverterType::ImageType InternalImage;
  typename InternalImage::Pointer img_ras = InternalImage::New();
  img_ras->Graft(image);

  m_Converter->SetVariable(varname, img_ras);
}

template <class TPixel, unsigned int VDim>
typename ConvertAPI<TPixel,VDim>::ImageType *
ConvertAPI<TPixel,VDim>
::GetImage(const char *varname)
{
  return m_Converter->GetVariable(varname);
}

template <class TPixel, unsigned int VDim>
bool 
ConvertAPI<TPixel,VDim>
::Execute(const char *cmdline, ...)
{
  char buffer[8192];

  va_list args;
  va_start (args, cmdline);
  vsprintf (buffer, cmdline, args);
  va_end (args);

  int argc = 0;
  char **argv = split_commandline(buffer, &argc);

  if(!argv)
    {
    m_Error = "Error parsing the command line expression";
    return false;
    }
  
  try 
    {
    int rc = m_Converter->ProcessCommandLine(argc, argv);
    if(rc)
      {
      m_Error = "Convert3D returned a non-zero return code. Check output for errors";
      return false;
      }
    else
      {
      return true;
      }
    }
  catch(std::exception &exc)
    {
    m_Error = exc.what();
    return false;
    }
  catch(...)
    {
    m_Error = "Unspecified failure. Check output for more information.";
    return false;
    }
}

template class ConvertAPI<double, 2>;
template class ConvertAPI<double, 3>;
template class ConvertAPI<double, 4>;
