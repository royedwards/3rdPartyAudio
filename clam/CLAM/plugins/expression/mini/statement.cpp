/*=============================================================================
    Copyright (c) 2001-2010 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#include "statement.hpp"

// This is not really called. Its only purpose is to
// instantiate the constructor of the grammar.
void instantiate_statement()
{
    typedef std::string::const_iterator iterator_type;
    functions_t functions;
    code_t code;
    statement<iterator_type> g(code, functions);
}
