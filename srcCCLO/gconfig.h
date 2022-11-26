/////////////////////////////////////////////////////////////////////////////
//////////          Evolutionary Algorithm (Configuration)        ///////////
//
// 27 Jun 2006, gkampolis
/////////////////////////////////////////////////////////////////////////////
//

#ifndef GCONFIG_H
#define GCONFIG_H

#include <iostream>

namespace NCONFIG
{

	
	class GConfig
	{
		protected:
			std::ostream& log;
		public:
			GConfig(std::ostream& _log);
			virtual ~GConfig(){};
		public:
			virtual bool ReadConfig(std::istream&)=0;
			virtual bool CheckConfig()=0;
	};



		
};
		
			

#endif
