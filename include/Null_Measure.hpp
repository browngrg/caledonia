#ifndef NULL_MEASURE_HPP
#define NULL_MEASURE_HPP


#include<string>
#include<stdio.h>


class Null_Measure
{
public:

   int    icall;         // number of calls
   bool   written;       // write since last sample?
   std::string header;   // Canned header information from other objects

public:

   Null_Measure();

   ~Null_Measure();

   // Add a string to the header information
   void add_header(std::string new_lines);

   // Add a sample state of walker to statistics (called by sampling object)
   template<typename MCWalker>
   void add_sample(MCWalker& walker, bool at_equilibrium=true);

   // Clear accummulated data and start over
   void clear();

   // Write current statistics (could be a check-point)
   void write();

   template<typename OBJ> void init(const OBJ& obj) {}

private:

};



Null_Measure::Null_Measure() 
{ 
   written=true; 
   icall = 0;
   // Any automatic set-up actions go here:
}

Null_Measure::~Null_Measure() 
{ 
   if( !written ) this->write(); 
   // Any clean-up actions go here:
}

void Null_Measure::add_header(std::string new_lines)
{
   header = header + new_lines;
}

template<typename MCWalker>
void Null_Measure::add_sample(MCWalker& walker, bool at_equilibrium)
{
   icall++;
   written = false;
   // Actions that happen regardless of equilibrium:
   // Equilibrium Measurements Come Here:
}

void Null_Measure::clear()
{
   written = false;
   icall = 0;
   // Actions to clear accumulated statistics go here:
}

void Null_Measure::write()
{
   written = true;
}

#endif
