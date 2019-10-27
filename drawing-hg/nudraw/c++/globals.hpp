#ifndef GLOBALS_HPP
#define GLOBALS_HPP

namespace nudraw
{
    namespace constants
    {
        // possible base types (internal representation)
        enum base_type {STRAND_START, STRAND_END, PAIR_LEFT, PAIR_RIGHT, UNPAIRED};
        
        // symbols for the various base types (strings to look for in parser)
        // XXX: auto-casts chars into ints
        static const char NICK_SYMBOL = '+';
        static const char LEFT_PAIR_SYMBOL = '(';
        static const char RIGHT_PAIR_SYMBOL = ')';
        static const char UNPAIRED_SYMBOL = '.';
    }
}

#endif        //  #ifndef GLOBALS_HPP

