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

#ifndef __Documentation_h__
#define __Documentation_h__

#include <string>
#include <set>
#include <vector>

/**
 * Documentation parsing system - takes our markdown format file and uses it
 * to generate help output
 */
class Documentation
{
public:

  struct CommandDoc
    {
    std::string Title;
    std::vector<std::string> Aliases;
    std::string ShortDesc;
    std::string LongDesc;
    };

  struct Category
    {
    std::string Title;
    std::vector<CommandDoc> Commands;
    Category(const std::string name) : Title(name) {}
    };

  /** Constructor - takes a C string with markdown text */
  Documentation(unsigned char* rawdoc);

  void PrintCommandListing(std::ostream &out);
  bool PrintCommandHelp(std::ostream &out, const std::string &command);
  void PrintManual(std::ostream &out);
  
  const std::set<std::string> &GetAllCommands() const
    { return m_AllCommands; }

  const std::vector<Category> &GetCategories() const
    { return m_Categories; }

  // Text processing routines
  static std::string &ltrim(std::string &s);
  static std::string &rtrim(std::string &s);
  static std::string &trim(std::string &s);
  static std::string tolower(const std::string &s);

private:
  // Complete manual text
  std::string m_Text;

  // Headings for commands and categories
  std::string m_CategoryHeading, m_CommandHeading;

  // Grouping of commands
  std::vector<Category> m_Categories;

  // Complete listing of commands and aliases
  std::set<std::string> m_AllCommands;
};




#endif

