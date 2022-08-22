/** \brief Basic includes and defining Log
 *
 *  Includes basic libraries.
 *  Sets up a Logarithm function that has size_t as its input and output type
 */
/*
   Copyright 2022 UChicago Argonne, LLC and contributors

   Licensed under the Apache License, Version 2.0 (the "License"); 
   you may not use this file except in compliance with the License. 
    
 
   Unless required by applicable law or agreed to in writing, software 
   distributed under the License is distributed on an "AS IS" BASIS, 
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
   See the License for the specific language governing permissions and 
   limitations under the License.
*/

#ifndef BITTREE_PRELUDE_H__
#define BITTREE_PRELUDE_H__

#include <cstddef>
#include <climits>
#include <stdexcept>
#include <memory>
#include <vector>
#include <string>

namespace bittree {

  /** Arbitrary base logarithm that takes size_t as input and output */ 
  template<size_t b, size_t x>
  struct Log {
               /** Recursive function to compute Log */
               enum : unsigned { val = 1u + Log<b,x/b>::val };
             };
  template<size_t b> struct Log<b,1> { enum { val = 0u }; };
  template<size_t b> struct Log<b,0> { enum { val = 0u }; };

}
#endif
