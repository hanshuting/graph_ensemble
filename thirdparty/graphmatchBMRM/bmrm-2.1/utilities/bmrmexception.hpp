#ifndef _BMRMEXCEPTION_HPP_
#define _BMRMEXCEPTION_HPP_

#include <exception>
#include <string>

using namespace std;

class CBMRMException : public exception {
   public:
      CBMRMException(const string& theMessage, const string& theThrower);
      virtual ~CBMRMException() throw();
      virtual const string& Report();

   protected:
      string message;
      string thrower;      
};


#endif
