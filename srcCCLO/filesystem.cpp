#include "filesystem.h"
#define _USEDIRECT3C
#ifdef WIN32
	#include <direct.h>
	#include "windows.h"
	#include <conio.h>
#else
	#ifdef _USEDIRECT3C
		#include <sys/types.h>
		#include <dirent.h>
	#endif
	#include <unistd.h>
	#include <sys/stat.h>
//	#include <curses.h>
#endif

using namespace std;

std::string GetPassword()
{
	std::string pass="";
	std::cin>>pass;
	return pass;
}


bool changedir(const std::string& dir)
{
#ifdef WIN32
	if (_chdir(dir.c_str())==0) return true;
	return false;
#else
	if (chdir(dir.c_str())==0) return true;
	return false;
#endif
}


bool makedir(const std::string& dir)
{
#ifdef WIN32
	if (_mkdir(dir.c_str())==0) return true;
	return false;
#else
	if (mkdir(dir.c_str(),0755)==0) return true;
	return false;
#endif
}


bool getworkdir(std::string& dir)
{
	dir="";
	const int maxlen=512;
	char buf[maxlen];	//A little big
#ifdef WIN32
	if (_getcwd(buf,maxlen)==NULL) return false;
#else
	if (getcwd(buf,maxlen)==NULL) return false;
#endif
	dir=buf;
	return true;
}


bool istream_operational(const std::istream& in)
{
	std::ios_base::iostate s=in.rdstate();
	if((s & std::ios::badbit)  ||
	   (s & std::ios::failbit) ||
	   (s & std::ios::eofbit) ) return false;
	return true;
}

void istream_next_char(std::istream& in, const char t)
{
	char p='\0';
	while ( istream_operational(in)&&p!=t)
		p=in.get();
}

void istream_nl(std::istream& in)
{
	istream_next_char(in,'\n');
}

void istream_get_tokens(std::istream& in, std::vector<std::string>& tokens)
{
	tokens.resize(0);
	std::string t="";
	char p='\0';
	while (istream_operational(in)&&p!='\n')
	{
		p=in.get();
		if (p!=' '&&p!='\n')t+=p;
		else
		{
			tokens.push_back(t);
			t="";
		}
	}
	if (t.size()>0) tokens.push_back(t);
}

void skipLines(std::istream& in, const int nLines)
{
	for (int i=0;i<nLines;i++) istream_nl(in);
}

int SystemCall(const std::string script)
{
	int sys_ret=0;

#ifdef WIN32
	PROCESS_INFORMATION proc_info;
	STARTUPINFO startup;
	startup.cb = sizeof(startup);
	startup.lpReserved = NULL;
	startup.lpDesktop = NULL;
	startup.lpTitle = NULL;
	startup.dwFlags = STARTF_USESHOWWINDOW;
	startup.wShowWindow = SW_HIDE;
	startup.cbReserved2 = 0;
	startup.lpReserved2 = NULL;

	//The following copy is done because the second argument of
	//CreateProcess is char* and not const char*
	char* cscript=new char[script.length()+1];
	strncpy(cscript,script.c_str(),script.length());
	cscript[script.length()]=0;
//	cscript[script.copy(cscript,std::string::npos)]=0;
        BOOL successproc = CreateProcess(0,cscript,0,0,false,
		CREATE_DEFAULT_ERROR_MODE,0,0,&startup,&proc_info);
	if(successproc)
	{
		WaitForSingleObject(proc_info.hProcess,
			 INFINITE);
		CloseHandle(proc_info.hThread);
		CloseHandle(proc_info.hProcess);
		sys_ret = 0;
	}
	else
	{
		// sys_ret = GetLastError();
		sys_ret = 1;
		 LPVOID lpMsgBuf;
		FormatMessage(
			FORMAT_MESSAGE_ALLOCATE_BUFFER |
			FORMAT_MESSAGE_FROM_SYSTEM |
			FORMAT_MESSAGE_IGNORE_INSERTS,
			NULL,
			GetLastError(),
			MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
			// Default language
			(LPTSTR) &lpMsgBuf,
			0,
			NULL
		);
		std::string msg =
			std::string("Error executing ")+
			script+
			std::string(":\n") +
			std::string((char*)lpMsgBuf);
		LocalFree( lpMsgBuf );
	//	throw Default_exception(msg);
		// throw Ceval_exception("Error executing task.bat");
	}
	delete[] cscript;
#else
	sys_ret = system(script.c_str());
#endif

	return sys_ret;

}

std::string FileToString(const std::string name)
{
	std::ifstream in(name.c_str(),std::ios::in);
	if (!in) return "";
	return std::string(std::istreambuf_iterator<char>(in),
		        std::istreambuf_iterator<char>());
}

std::string StreamToString(std::istream& in)
{
	return std::string(std::istreambuf_iterator<char>(in),
		        std::istreambuf_iterator<char>());
}

#ifndef WIN32
bool deletedir(const char *path)
{
#ifndef _USEDIRECT3C
	std::string ss=std::string("rm -rf ")+std::string(path);
	system(ss.c_str());
	return true;	//Cannot check if success
#else
	struct stat sb1;
     	if (-1 == stat(path, &sb1)) return false;
	if (!(S_IFDIR & sb1.st_mode))
	{
		// try !S_ISDIR(sb.st_mode)
		return false;
      	}
	//
	DIR *dir;
	if (NULL == (dir = opendir(path))) return false;
	if (-1 == chdir(path)) return false;
	struct dirent *dirent;
	struct stat sb;
	while (NULL != (dirent = readdir(dir)))
	{
		if (-1 == stat(dirent->d_name, &sb)) return false;
		if (strcmp(dirent->d_name, ".") &&
		    strcmp(dirent->d_name, ".."))
	       	{
			if (S_ISDIR(sb.st_mode))
			{
				if (!deletedir(dirent->d_name)) return false;
			}
			else if (-1 == unlink(dirent->d_name)) return false;
		}
	}
	if (-1 == chdir("..")) return false;
	if (-1 == closedir(dir)) return false;
	if (-1 == rmdir(path)) return false;
	return true;
#endif
}
#endif

