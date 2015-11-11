#ifndef PROGRAM_OPTIONS_H
#define PROGRAM_OPTIONS_H

// Class for setting parameters with default values, input file, or command line argument.
// Sept 16, 2014
// Greg Brown (browngrg@comacast.net,gbrown@fsu.edu)

// TODO: 8/6/2014  Need to store unrecognized options.
//       When get a new add_option, check to see if already read. If so, then set the value!!.
//       That is what map_unknown should be used for
//
//       2/26/2015 Added system information like start time for automatic reporting

#include<unistd.h>
#include<iostream>
#include<fstream>
#include<string>
#include<string.h>
#include<map>
#include<vector>

#include <cstring>
#include <cstdlib>
#include <getopt.h>
#include <stdio.h>


#define PROGRAMOPTIONS_USE_TIME
#ifdef PROGRAMOPTIONS_USE_TIME
#include<ctime>
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

struct POREC
{
public:
   std::string name;
   std::string desc;
   char short_opt;
   int  indx_global,indx_type;
   enum TYPE { t_bool, t_int, t_long, t_float, t_double, t_string } t_type;
public:
   POREC(std::string n, std::string d, char c, TYPE t, int g, int i) 
   { 
      name=n; 
      desc=d; 
      short_opt=c; 
      t_type = t;
      indx_global = g;
      indx_type = i;
   }
   POREC(const POREC& r)
   {
      name = r.name;
      desc = r.desc;
      short_opt = r.short_opt;
      t_type = r.t_type;
      indx_global = r.indx_global;
      indx_type = r.indx_type;
   }
};  


// Read options from input file or command line
//
// Adds automation to application input without creating dependencies in code.
// Allows direct modification of variables and data members by associating "name" with pointer.
// Outputs a list of all parameters set when command line is pared, object is destructed, or write(fout) is called.
// This list can be edited, and then specified for input with the -i filename or --input=filename options
// The file is read when the option is encountered, subsequent command line optoins will override.
//
// Based on earlier versions of ProgramOptions and GenProgOptions
// from psimag toolset, altough none of the original code in included.
// The original concept is due to Thomas Schulthess.
//
class ProgramOptions
{
public:
   ProgramOptions(std::string name, std::string purpose="");
   ~ProgramOptions();
   void parse_command_line2(int argc, char* argv[]);
   void parse_command_line(int argc, char* argv[]);
   void read(std::istream& fin);
   void write();
   void write_help(std::ostream& fout);
   void mpi_synch();
public:
   std::vector<std::string> _argv;
   int process_id;
   int mpi_iproc;
   int mpi_nproc;
   bool mpi_enabled;
   std::string start_time_str;
   std::string writ_time_str;
#  ifdef PROGRAMOPTIONS_USE_TIME
   time_t start_time;
   time_t writ_time;
#  else
   int start_time;
   int writ_time;
#  endif
public:
   void add_option(std::string name, std::string descript, char schar, bool*   variable);
   void add_option(std::string name, std::string descript, char schar, int*    variable);
   void add_option(std::string name, std::string descript, char schar, long*   variable);
   void add_option(std::string name, std::string descript, char schar, float*  variable);
   void add_option(std::string name, std::string descript, char schar, double* variable);
   void add_option(std::string name, std::string descript, char schar, char*   variable);
   void get_value(std::string name, bool& value) const;
   void get_value(std::string name, int& value) const;
   void get_value(std::string name, long& value) const;
   void get_value(std::string name, float& value) const;
   void get_value(std::string name, double& value) const;
   void get_value(std::string name, std::string& value) const;
   template<typename T> T get_value(std::string name) const;
   void set_value(std::string name, bool value);
   void set_value(std::string name, int  value);
   void set_value(std::string name, long value);
   void set_value(std::string name, float  value);
   void set_value(std::string name, double value);
   void set_value(std::string name, std::string value);
private:
   typedef std::map<std::string,int>::const_iterator POPTR;
   std::map<std::string,int>   map_option;
   int                         char_option[128];
   std::vector<POREC>          option;
   std::vector<bool*>          _var_bool;
   std::vector<int*>           _var_int;
   std::vector<long*>          _var_long;
   std::vector<float*>         _var_float;
   std::vector<double*>        _var_double;
   std::vector<char*>          _var_string;
private:
   std::map<std::string,std::string> map_unknown;
private:
   std::string _prog,_purp;
   char        _fname_in[128],_fname_out[128];
   bool        _help,_verbose;
};


ProgramOptions::ProgramOptions(std::string name, std::string purpose)
{
   _prog = name;
   _purp = purpose;
   // clear flag map
   for(int i=0; i<128; i++) char_option[i] = -1;
   // default options
   add_option("help","enable this output",'h',&_help);
   add_option("verbose","enable verbose output",'v',&_verbose);
   add_option("input","filename for options input",'i',_fname_in);
   add_option("output","filename for options output",'o',_fname_out);
   // default values
   _help      = false;
   _verbose   = false;
   sprintf(_fname_in,"");
   sprintf(_fname_out,"ProgramOptionsOut.txt"); 
   mpi_enabled = false;
   mpi_iproc = 0;
   mpi_nproc = 1;
   start_time = 0;
   start_time_str = "";
   writ_time = 0;
   writ_time_str = "unfinished";
   // get run information
   process_id = getpid();
#  ifdef PROGRAMOPTIONS_USE_TIME
   time(&start_time);
   struct tm * timeinfo = localtime(&start_time);
   start_time_str = asctime(timeinfo);
   start_time_str = start_time_str.substr(0,start_time_str.size()-1);
#  endif
#  ifdef USE_MPI
   mpi_enabled = true;
   MPI_Comm_rank(MPI_COMM_WORLD,&mpi_iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&mpi_nproc);
#  endif
}


ProgramOptions::~ProgramOptions()
{
}


void ProgramOptions::set_value(std::string name, bool value)
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second;
   char buffer[100];
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   *(_var_bool  [option[iopt].indx_type]) = static_cast<bool>(value); break;
   case POREC::t_int:    *(_var_int   [option[iopt].indx_type]) = static_cast<int>(value); break;
   case POREC::t_long:   *(_var_long  [option[iopt].indx_type]) = static_cast<long>(value); break;
   case POREC::t_float:  *(_var_float [option[iopt].indx_type]) = static_cast<float>(value); break;
   case POREC::t_double: *(_var_double[option[iopt].indx_type]) = static_cast<double>(value); break;
   case POREC::t_string: 
      sprintf(buffer,"%d",static_cast<int>(value)); 
      strcpy(_var_string[option[iopt].indx_type],buffer);
      break;
   }
}


void ProgramOptions::set_value(std::string name, int value)
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second;
   char buffer[100];
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   *(_var_bool  [option[iopt].indx_type]) = static_cast<bool>(value); break;
   case POREC::t_int:    *(_var_int   [option[iopt].indx_type]) = static_cast<int>(value); break;
   case POREC::t_long:   *(_var_long  [option[iopt].indx_type]) = static_cast<long>(value); break;
   case POREC::t_float:  *(_var_float [option[iopt].indx_type]) = static_cast<float>(value); break;
   case POREC::t_double: *(_var_double[option[iopt].indx_type]) = static_cast<double>(value); break;
   case POREC::t_string:
      sprintf(buffer,"%d",value); 
      strcpy(_var_string[option[iopt].indx_type],buffer);
      break;
   }
}


void ProgramOptions::set_value(std::string name, long value)
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second;
   char buffer[100];
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   *(_var_bool  [option[iopt].indx_type]) = static_cast<bool>(value); break;
   case POREC::t_int:    *(_var_int   [option[iopt].indx_type]) = static_cast<int>(value); break;
   case POREC::t_long:   *(_var_long  [option[iopt].indx_type]) = static_cast<long>(value); break;
   case POREC::t_float:  *(_var_float [option[iopt].indx_type]) = static_cast<float>(value); break;
   case POREC::t_double: *(_var_double[option[iopt].indx_type]) = static_cast<double>(value); break;
   case POREC::t_string:
      sprintf(buffer,"%ld",value); 
      strcpy(_var_string[option[iopt].indx_type],buffer);
      break;
   }
}


void ProgramOptions::set_value(std::string name, float value)
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second;
   char buffer[100];
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   *(_var_bool  [option[iopt].indx_type]) = static_cast<bool>(value); break;
   case POREC::t_int:    *(_var_int   [option[iopt].indx_type]) = static_cast<int>(value); break;
   case POREC::t_long:   *(_var_long  [option[iopt].indx_type]) = static_cast<long>(value); break;
   case POREC::t_float:  *(_var_float [option[iopt].indx_type]) = static_cast<float>(value); break;
   case POREC::t_double: *(_var_double[option[iopt].indx_type]) = static_cast<double>(value); break;
   case POREC::t_string:
      sprintf(buffer,"%f",value); 
      strcpy(_var_string[option[iopt].indx_type],buffer);
      break;
   }
}


void ProgramOptions::set_value(std::string name, double value)
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second;
   char buffer[100];
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   *(_var_bool  [option[iopt].indx_type]) = static_cast<bool>(value); break;
   case POREC::t_int:    *(_var_int   [option[iopt].indx_type]) = static_cast<int>(value); break;
   case POREC::t_long:   *(_var_long  [option[iopt].indx_type]) = static_cast<long>(value); break;
   case POREC::t_float:  *(_var_float [option[iopt].indx_type]) = static_cast<float>(value); break;
   case POREC::t_double: *(_var_double[option[iopt].indx_type]) = static_cast<double>(value); break;
   case POREC::t_string:
      sprintf(buffer,"%lf",value); 
      strcpy(_var_string[option[iopt].indx_type],buffer);
      break;
   }
}


void ProgramOptions::set_value(std::string name, std::string value)
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second;
   int itmp;
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   sscanf( value.c_str(), "%d", &itmp); *(_var_bool[ option[iopt].indx_type ]) = itmp; break;
   case POREC::t_int:    sscanf( value.c_str(), "%d",  _var_int   [ option[iopt].indx_type ] ); break;
   case POREC::t_long:   sscanf( value.c_str(), "%ld", _var_long  [ option[iopt].indx_type ] ); break;
   case POREC::t_float:  sscanf( value.c_str(), "%f",  _var_float [ option[iopt].indx_type ] ); break;
   case POREC::t_double: sscanf( value.c_str(), "%lf", _var_double[ option[iopt].indx_type ] ); break;
   case POREC::t_string: strcpy(_var_string[option[iopt].indx_type],value.c_str()); break;
   }
}


void ProgramOptions::get_value(std::string name, bool& value) const
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second; 
   int iflag;
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   value = static_cast<bool>( *(_var_bool  [option[iopt].indx_type]) ); break;
   case POREC::t_int:    value = static_cast<bool>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_long:   value = static_cast<bool>( *(_var_long  [option[iopt].indx_type]) ); break;
   case POREC::t_float:  value = static_cast<bool>( *(_var_float [option[iopt].indx_type]) ); break;
   case POREC::t_double: value = static_cast<bool>( *(_var_double[option[iopt].indx_type]) ); break;
   case POREC::t_string: value = sscanf(_var_string[option[iopt].indx_type],"%d",&iflag); value=static_cast<bool>(iflag); break;
   }
}


void ProgramOptions::get_value(std::string name, int& value) const
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second; 
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   value = static_cast<int>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_int:    value = static_cast<int>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_long:   value = static_cast<int>( *(_var_long  [option[iopt].indx_type]) ); break;
   case POREC::t_float:  value = static_cast<int>( *(_var_float [option[iopt].indx_type]) ); break;
   case POREC::t_double: value = static_cast<int>( *(_var_double[option[iopt].indx_type]) ); break;
   case POREC::t_string: value = sscanf(_var_string[option[iopt].indx_type],"%d",&value); break;
   }
}


void ProgramOptions::get_value(std::string name, long& value) const
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second; 
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   value = static_cast<long>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_int:    value = static_cast<long>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_long:   value = static_cast<long>( *(_var_long  [option[iopt].indx_type]) ); break;
   case POREC::t_float:  value = static_cast<long>( *(_var_float [option[iopt].indx_type]) ); break;
   case POREC::t_double: value = static_cast<long>( *(_var_double[option[iopt].indx_type]) ); break;
   case POREC::t_string: value = sscanf(_var_string[option[iopt].indx_type],"%ld",&value); break;
   }
}


void ProgramOptions::get_value(std::string name, float& value) const
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second; 
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   value = static_cast<float>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_int:    value = static_cast<float>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_long:   value = static_cast<float>( *(_var_long  [option[iopt].indx_type]) ); break;
   case POREC::t_float:  value = static_cast<float>( *(_var_float [option[iopt].indx_type]) ); break;
   case POREC::t_double: value = static_cast<float>( *(_var_double[option[iopt].indx_type]) ); break;
   case POREC::t_string: value = sscanf(_var_string[option[iopt].indx_type],"%f",&value); break;
   }
}


void ProgramOptions::get_value(std::string name, double& value) const
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second; 
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   value = static_cast<double>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_int:    value = static_cast<double>( *(_var_int   [option[iopt].indx_type]) ); break;
   case POREC::t_long:   value = static_cast<double>( *(_var_long  [option[iopt].indx_type]) ); break;
   case POREC::t_float:  value = static_cast<double>( *(_var_float [option[iopt].indx_type]) ); break;
   case POREC::t_double: value = static_cast<double>( *(_var_double[option[iopt].indx_type]) ); break;
   case POREC::t_string: value = sscanf(_var_string[option[iopt].indx_type],"%lf",&value); break;
   }
}


void ProgramOptions::get_value(std::string name, std::string& value) const
{
   POPTR optr = map_option.find(name);
   if( optr==map_option.end() ) return;
   int iopt = optr->second;
   char buffer[100];
   switch(option[iopt].t_type)
   {
   case POREC::t_bool:   sprintf(buffer,"%d",static_cast<int>(*(_var_bool[option[iopt].indx_type]))); break;
   case POREC::t_int:    sprintf(buffer,"%d",*(_var_int[option[iopt].indx_type])); break;
   case POREC::t_long:   sprintf(buffer,"%ld",*(_var_long[option[iopt].indx_type])); break;
   case POREC::t_float:  sprintf(buffer,"%f",*(_var_float[option[iopt].indx_type])); break;
   case POREC::t_double: sprintf(buffer,"%lf",*(_var_double[option[iopt].indx_type])); break;
   case POREC::t_string: strcpy(buffer,_var_string[option[iopt].indx_type]); break;
   }
   value = std::string(buffer);
}


template<typename T>
T ProgramOptions::get_value(std::string name) const { T t; get_value(name,t); return t; }


void ProgramOptions::add_option(std::string name, std::string descript, char schar, bool*   variable)
{
   POPTR ptr = map_option.find(name);
   if( ptr!=map_option.end() ) return;
   POREC po(name,descript,schar,POREC::t_bool,option.size(),_var_bool.size());
   map_option.insert( std::make_pair(name,po.indx_global) );
   option.push_back(po);
   _var_bool.push_back(variable);
   char_option[schar] = po.indx_global;
}


void ProgramOptions::add_option(std::string name, std::string descript, char schar, int*    variable)
{
   POPTR ptr = map_option.find(name);
   if( ptr!=map_option.end() ) return;
   POREC po(name,descript,schar,POREC::t_int,option.size(),_var_int.size());
   map_option.insert( std::make_pair(name,po.indx_global) );
   option.push_back(po);
   _var_int.push_back(variable);
   char_option[schar] = po.indx_global;
}


void ProgramOptions::add_option(std::string name, std::string descript, char schar, long*    variable)
{
   POPTR ptr = map_option.find(name);
   if( ptr!=map_option.end() ) return;
   POREC po(name,descript,schar,POREC::t_long,option.size(),_var_long.size());
   map_option.insert( std::make_pair(name,po.indx_global) );
   option.push_back(po);
   _var_long.push_back(variable);
   char_option[schar] = po.indx_global;
}


void ProgramOptions::add_option(std::string name, std::string descript, char schar, float*  variable)
{
   POPTR ptr = map_option.find(name);
   if( ptr!=map_option.end() ) return;
   POREC po(name,descript,schar,POREC::t_float,option.size(),_var_float.size());
   map_option.insert( std::make_pair(name,po.indx_global) );
   option.push_back(po);
   _var_float.push_back(variable);
   char_option[schar] = po.indx_global;
}


void ProgramOptions::add_option(std::string name, std::string descript, char schar, double* variable)
{
   POPTR ptr = map_option.find(name);
   if( ptr!=map_option.end() ) return;
   POREC po(name,descript,schar,POREC::t_double,option.size(),_var_double.size());
   map_option.insert( std::make_pair(name,po.indx_global) );
   option.push_back(po);
   _var_double.push_back(variable);
   char_option[schar] = po.indx_global;
}


void ProgramOptions::add_option(std::string name, std::string descript, char schar, char* variable)
{
   POPTR ptr = map_option.find(name);
   if( ptr!=map_option.end() ) return;
   POREC po(name,descript,schar,POREC::t_string,option.size(),_var_string.size());
   map_option.insert( std::make_pair(name,po.indx_global) );
   option.push_back(po);
   _var_string.push_back(variable);
   char_option[schar] = po.indx_global;
}


void ProgramOptions::read(std::istream& fin)
{
   if( !fin ) return;
   while( !fin.eof() )
   {
      std::string linebuff;
      std::getline(fin,linebuff);
      if( linebuff.size()>0 && linebuff[0]!='#' )
      {
         size_t iptr = linebuff.find_first_of(":=");
         if( iptr != std::string::npos )
         {
            size_t jptr;
            std::string name = linebuff.substr(0,iptr);
            std::string val  = linebuff.substr(iptr+1);
            jptr = name.find_first_not_of(" \t");
            if(jptr!=std::string::npos) name = name.substr(jptr);
            jptr = name.find_last_not_of(" \t");
            if(jptr!=std::string::npos) name = name.substr(0,jptr+1);
            jptr = val.find_first_not_of(" \t");
            if(jptr!=std::string::npos) val = val.substr(jptr);
            jptr = val.find_last_not_of(" \t");
            if(jptr!=std::string::npos) val = val.substr(0,jptr+1);
            set_value(name,val);
         }
      }
   }
}


void ProgramOptions::write()
{
   if( mpi_iproc>0 ) return;
   std::string ofname = get_value<std::string>("output");
   std::ofstream fout(ofname.c_str());
   double secs = 0;
   double psecs = 0;
#  ifdef PROGRAMOPTIONS_USE_TIME
   time(&writ_time);
   struct tm * timeinfo = localtime(&writ_time);
   writ_time_str = asctime(timeinfo);
   writ_time_str = writ_time_str.substr(0,writ_time_str.size()-1);
   secs = difftime(writ_time,start_time);
   clock_t tics = clock();
   psecs = static_cast<double>(tics)/static_cast<double>(CLOCKS_PER_SEC);
#  endif
   int resident = 0;
#  if 1
   int tSize = 0, share = 0;                          // In kilobytes
   std::ifstream buffer("/proc/self/statm");
   buffer >> tSize >> resident >> share;
   buffer.close();
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   resident = resident * page_size_kb;
#   endif
   fout << "# " <<  _prog << ": " << _purp << std::endl;
   fout << "# process id: " << process_id << std::endl;
   fout << "# mpi enabled: " << mpi_enabled << std::endl;
   fout << "# mpi nproc: " << mpi_nproc << std::endl;
   fout << "# start time: " << start_time_str << std::endl;
   fout << "# write time: " << writ_time_str << std::endl;
   fout << "# wall time: " << secs << " seconds" << std::endl;
   fout << "# process time: " << psecs << " seconds" << std::endl;
   fout << "# resident memory = " << resident << " KBytes" << std::endl;
   share = share * page_size_kb;
   for(int i=0; i<option.size(); i++)
   {
      fout << option[i].name << "= ";
      switch(option[i].t_type)
      {
         case POREC::t_bool:   fout << *(_var_bool  [ option[i].indx_type ]); break;
         case POREC::t_int:    fout << *(_var_int   [ option[i].indx_type ]); break;
         case POREC::t_long:   fout << *(_var_long  [ option[i].indx_type ]); break;
         case POREC::t_float:  fout << *(_var_float [ option[i].indx_type ]); break;
         case POREC::t_double: fout << *(_var_double[ option[i].indx_type ]); break;
         case POREC::t_string: fout <<   _var_string[ option[i].indx_type ] ; break;
      }
      fout << std::endl;
   }
}


void ProgramOptions::parse_command_line2(int argc, char* argv[])
{
   // look for input file
   for(int iarg=1; iarg<argc; iarg++)
   {
      std::string arg(argv[iarg]);
      size_t strp = arg.find("--input=");
      if( strp!=std::string::npos )
      {
         strp += 7;
         strcpy(_fname_in,argv[iarg]+strp);
      }
      if( strcmp("-i",argv[iarg])==0 )
      {
         strcpy(_fname_in,argv[iarg+1]);
      } 
   }
   std::string fname = get_value<std::string>("input");
   std::ifstream fin(fname.c_str());
   read(fin);
   // parse command line
   for(int iarg=1; iarg<argc; iarg++)
   {
      std::string arg(argv[iarg]);
      if( arg.size()>2 && arg[0]=='-' )
      {
         if( arg[1]=='-' )
         {
            size_t strp = arg.find("=");
            std::string name = arg.substr(2,strp-2);
            set_value(name,arg.substr(strp+1));
         }
         else
         {
            int iopt = char_option[arg[1]];
            if( iopt=='v' ) set_value(option['v'].name,true);
            else if( iopt>=0 && iopt<128 )
            set_value(option[iopt].name,argv[iarg+1]);
         }
      }
   }
}


void ProgramOptions::parse_command_line(int argc, char* argv[])
{
   if( mpi_iproc==0 )
   {
   // Make a copy of the command line, and of unknown parameters
   _argv.resize(argc);
   for(int iarg=0; iarg<argc; iarg++) argv[iarg] = argv[iarg];
   map_unknown.clear();
   // Actual parsing uses this documentation:
   // name, has_arg, flag, val
   // const char *name
   //     This field is the name of the option. It is a string.
   // int has_arg
   //     This field says whether the option takes an argument. 
   //     It is an integer, and there are three legitimate values: 
   //     no_argument, required_argument and optional_argument.
   // int *flag
   // int val
   //     These fields control how to report or act on the option when it occurs.
   //     If flag is a null pointer, then the val is a value which identifies this option. 
   //     Often these values are chosen to uniquely identify particular long options.
   //     If flag is not a null pointer, it should be the address of an int variable 
   //     which is the flag for this option. The value in val is the value to store in 
   //     the flag to indicate that the option was seen. 
   int nopt = option.size();
   std::vector<struct option> long_options( nopt+1 );
   for(int i=0; i<nopt; i++)
   {
      long_options[i].name = option[i].name.c_str();
      if( option[i].t_type==POREC::t_bool )
         long_options[i].has_arg = no_argument;
      else
         long_options[i].has_arg = true;
      long_options[i].flag = 0;
      long_options[i].val = i + 128;
   }
   long_options[nopt].name    = 0;   // end of options record
   long_options[nopt].has_arg = 0;
   long_options[nopt].flag    = 0;
   long_options[nopt].val     = 0;
   // short options
   std::string short_options;
   for(int i=0; i<nopt; i++)
   {
      if( option[i].short_opt!=' ' )
      {
         short_options += option[i].short_opt;
         if( option[i].t_type!=POREC::t_bool )
            short_options += ':';
      }
   } 
   int input_opt = map_option["input"];   // Read from input file first
   while(true)
   {
      // Parse the command line for options
      int option_index = 0;
      int c = getopt_long(argc,argv,short_options.c_str(),&(long_options[0]),&option_index);
      if(c == -1) break;            // Detect the end of the options.
      int i = ((c>=0)&&(c<128))? char_option[c] : c - 128;
      switch( option[i].t_type )
      {
      case POREC::t_bool:   *(_var_bool[ option[i].indx_type ]) = true;                  break;
      case POREC::t_int:    sscanf( optarg, "%d",  _var_int   [ option[i].indx_type ] ); break;
      case POREC::t_long:   sscanf( optarg, "%ld", _var_long  [ option[i].indx_type ] ); break;
      case POREC::t_float:  sscanf( optarg, "%f",  _var_float [ option[i].indx_type ] ); break;
      case POREC::t_double: sscanf( optarg, "%lf", _var_double[ option[i].indx_type ] ); break;
      case POREC::t_string: sscanf( optarg, "%s",  _var_string[ option[i].indx_type ] ); break;
      }
      if( i==input_opt )
      {
         std::string fname = get_value<std::string>("input");
         std::ifstream fin(fname.c_str());
         read(fin);
      }
   }
   if( get_value<bool>("help")==true )
   {
      write_help(std::cout);
   }
   else
   {
      write();  
   }
   }
   mpi_synch();
   write();
}


void ProgramOptions::write_help(std::ostream& fout)
{
   fout << _prog << ": " << _purp << std::endl;
   for(int i=0; i<option.size(); i++)
   {
      if( option[i].short_opt!=' ' ) 
         fout << " -" << option[i].short_opt;
      else
         fout << "   ";
      fout << "  --" << option[i].name;
      switch(option[i].t_type)
      {
      case POREC::t_bool:   fout << "        "; break;
      case POREC::t_int:    fout << " int    "; break;
      case POREC::t_long:   fout << " long   "; break;
      case POREC::t_float:  fout << " float  "; break;
      case POREC::t_double: fout << " double "; break;
      case POREC::t_string: fout << " string "; break;
      }
      fout << "  " << option[i].desc << std::endl;
   }
}

void ProgramOptions::mpi_synch()
{
#  ifdef USE_MPI
   // Need to transfer the values in each of these
   // std::vector<bool*>          _var_bool;
   // std::vector<int*>           _var_int;
   // std::vector<long*>          _var_long;
   // std::vector<float*>         _var_float;
   // std::vector<double*>        _var_double;
   // std::vector<char*>          _var_string;
   int nstr = _var_string.size();
   for(int i=0; i<_var_string.size(); i++) nstr += std::strlen(_var_string[i]);
   {
      std::vector<long> buffer( _var_bool.size() + _var_int.size() + _var_long.size() + 1 );
      int ib = 0;
      for(int i=0; i<_var_bool.size(); i++) buffer[ib++] = *(_var_bool[i]);
      for(int i=0; i<_var_int.size();  i++) buffer[ib++] = *(_var_int[i]);
      for(int i=0; i<_var_long.size(); i++) buffer[ib++] = *(_var_long[i]);
      buffer[ib++] = nstr;
      int nbuff = buffer.size();
      MPI_Bcast(&(buffer[0]),nbuff,MPI_LONG,0,MPI_COMM_WORLD);
      ib = 0;
      for(int i=0; i<_var_bool.size(); i++) *(_var_bool[i]) = buffer[ib++];
      for(int i=0; i<_var_int.size();  i++) *(_var_int[i] ) = buffer[ib++];
      for(int i=0; i<_var_long.size(); i++) *(_var_long[i]) = buffer[ib++];
      nstr = buffer[ib++];
   }
   {
      std::vector<double> buffer( _var_float.size() + _var_double.size() );
      int ib = 0;
      for(int i=0; i<_var_float.size(); i++)  buffer[ib++] = *(_var_float[i]);
      for(int i=0; i<_var_double.size(); i++) buffer[ib++] = *(_var_double[i]);
      int nbuff = buffer.size();
      MPI_Bcast(&(buffer[0]),nbuff,MPI_DOUBLE,0,MPI_COMM_WORLD);
      ib = 0;
      for(int i=0; i<_var_float.size(); i++)  *(_var_float[i])  = buffer[ib++];
      for(int i=0; i<_var_double.size(); i++) *(_var_double[i]) = buffer[ib++];
   }
   {
      // nstr was synchronized with other integers
      std::vector<char> buffer(nstr);
      int ib = 0;
      for(int i=0; i<_var_string.size(); i++)
      {
         std::strcpy(&(buffer[ib]),_var_string[i]);
         ib += std::strlen(_var_string[i]) + 1;
      }
      MPI_Bcast(&(buffer[0]),buffer.size(),MPI_CHAR,0,MPI_COMM_WORLD);
      ib = 0;
      for(int i=0; i<_var_string.size(); i++)
      {
         std::strcpy(_var_string[i],&(buffer[ib]));
         ib += std::strlen(_var_string[i]) + 1;
      }
   }
#  endif
}

#endif   // PROGRAM_OPTIONS_H 
