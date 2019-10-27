#include <iostream>

#include "base.hpp"
#include "strand.hpp"

#include <boost/program_options.hpp>
#include <iterator>

using namespace std;

int main(int argc, char* argv[])
{
    try
    {
        //boost::program_options::options_description desc("Allowed options");
        printf("In a try block \n");
//        desc.add_options()("help", "produce help message")
        
    }
    
    catch(exception& e)
    {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    
    catch(...)
    {
        cerr << "An uncaught exception was thrown\n";
        return 1;
    }
    
    
    nudraw::base *mybase = new nudraw::base(nudraw::constants::PAIR_LEFT, 0);
    
    cout << "type: " << mybase->get_type() << endl;
    cout << "index: " << mybase->get_index() << endl;
    //std::cout << "x pos: " << mybase->get_position().get_x() << std::endl;
    
    delete mybase;
    
    return 0;
}
