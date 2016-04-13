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

#include "Documentation.h"
#include <iomanip>
#include <sstream>
#include <itksys/RegularExpression.hxx>
#include <algorithm>
#include <functional>
#include <cctype>

// Helper function: trim whitespace
// from: http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
//
// trim from start
std::string &
Documentation::ltrim(std::string &s) 
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
std::string &
Documentation::rtrim(std::string &s) 
{
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
std::string &
Documentation::trim(std::string &s) 
{
  return ltrim(rtrim(s));
}

std::string 
Documentation::tolower(const std::string &src)
{
  std::string s = src;
  std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) ::tolower);
  return s;
}

Documentation::Documentation(unsigned char *rawdoc)
: m_Text((const char *) rawdoc)
{
  // Headings
  std::string relevant_category_heading = "### Commands:";
  m_CategoryHeading = "### ";
  m_CommandHeading = "#### ";

  // Parse the file
  std::istringstream iss(m_Text);

  // Current category
  int current_cat = -1;
  std::string current_command = "";

  // Read lines from the file
  std::string line;
  while(std::getline(iss, line)) 
    {
    if(line.find(m_CategoryHeading) == 0)
      {
      // This is a category - find out which one
      current_cat = -1;

      if(line.find(relevant_category_heading) == 0)
        {
        std::string title = line.substr(relevant_category_heading.length());
        m_Categories.push_back(Category(trim(title)));
        current_cat = ((int) m_Categories.size()) - 1;
        current_command = "";
        }
      }

    else if(line.find(m_CommandHeading) == 0 && current_cat >= 0)
      {
      // This is a parseable command. Parse out all of its information
      std::string rexp1 = m_CommandHeading + " *(.*): *(.*)$";
      itksys::RegularExpression re1(rexp1.c_str());
      if(re1.find(line))
        {
        CommandDoc cdoc;
        cdoc.Title = re1.match(1);
        cdoc.ShortDesc = re1.match(2);
        current_command = cdoc.Title;

        // Split the title into aliases
        std::istringstream isstmp(cdoc.Title);
        std::string cmdline;
        while(std::getline(isstmp, cmdline, ','))
          {
          std::string cmd = trim(cmdline);
          cdoc.Aliases.push_back(cmd);
          m_AllCommands.insert(cmd);
          }

        m_Categories[current_cat].Commands.push_back(cdoc);
        }
      else
        {
        current_command = "";
        }
      }

    else if(current_command.length() > 0 && current_cat >= 0)
      {
      // Line from a description
      m_Categories[current_cat].Commands.back().LongDesc += line;
      m_Categories[current_cat].Commands.back().LongDesc += "\n";
      }
    } 
}

void Documentation::PrintCommandListing(std::ostream &out)
{
  for(int i = 0; i < m_Categories.size(); i++)
    {
    out << m_Categories[i].Title << ": " << std::endl;
    for(int j = 0; j < m_Categories[i].Commands.size(); j++)
      {
      out << "    ";
      out << std::setw(32) << std::left;
      out << m_Categories[i].Commands[j].Title;
      out << ": ";
      out << m_Categories[i].Commands[j].ShortDesc;
      out << std::endl;
      }
    }
}

void Documentation::PrintManual(std::ostream &out)
{
  out << m_Text << std::endl;
}

bool Documentation::PrintCommandHelp(std::ostream &out, const std::string &command)
{
  // Create a search string
  if(command.length() == 0)
    return false;

  std::string req = command[0] == '-' ? command : std::string("-") + command;
  req = tolower(req);

  for(int i = 0; i < m_Categories.size(); i++)
    {
    for(int j = 0; j < m_Categories[i].Commands.size(); j++)
      {
      CommandDoc &cmd = m_Categories[i].Commands[j];
      for(int k = 0; k < cmd.Aliases.size(); k++)
        {
        if(tolower(cmd.Aliases[k]) == req)
          {
          out << std::setw(32) << std::left;
          out << m_Categories[i].Commands[j].Title;
          out << ": ";
          out << m_Categories[i].Commands[j].ShortDesc;
          out << std::endl;
          out << cmd.LongDesc;
          out << std::endl;
          return true;
          }
        }
      }
    }

  return false;
}
