/** \brief Basic includes and defining Log
 *
 *  Includes basic libraries.
 *  Sets up a Logarithm function that has size_t as its input and output type
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
