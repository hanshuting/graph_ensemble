#ifndef _COMMON_CPP_
#define _COMMON_CPP_ 

#include <string>
#include <vector>
#include <cctype>

/** Check if a string is blank.
 *
 *  \param line [read] The line to be checked if it is blank
 */
bool IsBlankLine(std::string &line)
{
        size_t n = line.size();
        for(size_t i=0; i<n; i++)
                if(!std::isspace(line[i]))
                        return false;
        return true;
}



void trim(std::string& str)
{  
        typedef std::string::size_type str_pos;
        str_pos pos = str.find_last_not_of(" \t\r");
        if(pos != std::string::npos) 
        {
                str.erase(pos + 1);
                pos = str.find_first_not_of(" \t\r");
                if(pos != std::string::npos) 
                        str.erase(0, pos);
        }
        else 
                str.erase(str.begin(), str.end());
}


void tokenize(const std::string& str, 
              std::vector<std::string>& tokens, 
              const std::string& delimiter = " ") {
  
  // svnvish: BUGBUG
  // delimiter is a character or a string?
  typedef std::string::size_type str_pos;
  
  // Skip delimiter at beginning.
  str_pos lastPos = str.find_first_not_of(delimiter, 0);
  // Find first "non-delimiter".
  str_pos pos = str.find_first_of(delimiter, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos){
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiter.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiter, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiter, lastPos);
  }
  return;
  
}

#endif
