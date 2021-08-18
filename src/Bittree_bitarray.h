#ifndef BITTREE_BITARRAY_H__
#define BITTREE_BITARRAY_H__

#include "Bittree_prelude.h"

namespace bittree {

  // Use unsigned int for word type
  typedef unsigned int WType;


  /** Stores, reads, and writes Bit Arrays.
   *    */
  class BitArray {
  /** W is the "word" type, looking at bittree_core its always instantiated to
   *  "unsigned". This is not incorrect but probably less efficient than using
   *  std::uintptr_t which will always be the CPUs native register bitwidth
   *  and is basically always 64-bits in HPC.

   *  BitArray works by packing the 1's and 0's into the binary representaiton of
   *  "words". A 32-bit integer can hold 32 bits. The least significant bit is
   *  in position zero. */
  
  public:
    /** Log_2 of memory allocated to objects of class W, in bits */
    static const unsigned logw = Log<2,sizeof(WType)*CHAR_BIT>::val; 
    /** bitwidth of WType */
    static const unsigned bitw = 1u << logw;
    static const WType one = WType(1);    /**< 1 cast as WType */
    static const WType ones = ~WType(0);  /**< Maximum length string of binary 1s cast as WType */
  private:
    unsigned len;       /**< Length of Bit Array. Access with length() */
    std::vector<WType> wbuf;          /**< Word buffer of type WType. Access with word_buf() */
  public:
    BitArray(unsigned len): len(len) {}
  public:
    static std::shared_ptr<BitArray > make(unsigned n);
  public:
    unsigned length() const { return len; }
    unsigned word_count() const { return (len+bitw-1u)>>logw; }
    std::vector<WType> word_buf() { return wbuf; }
    
    bool get(unsigned ix) const;     /**< get value of bit ix */
    bool set(unsigned ix, bool x);   /**< set value of bit ix */
    
    /** count 1's in interval [ix0,ix1) */
    unsigned count(unsigned ix0, unsigned ix1) const;
    /** count 1's in whole array */
    unsigned count() const;
    /** count 1's in either a or b */
    static unsigned count_xor(
      const std::shared_ptr<BitArray> a,
      const std::shared_ptr<BitArray> b,
      unsigned ix0, unsigned ix1
    );
    
    unsigned find(unsigned ix0, unsigned nth) const;
    
    void fill(bool x);                               /** Fill whole Bit Array */
    void fill(bool x, unsigned ix0, unsigned ix1);   /** Fill part of Bit Array */
    
    /** Class for reading off the bit array.
     *  Takes a virtual BitArray of type W and reads off values.
     *  Note: reading off the end is safe and returns 0 bits   */
    class Reader {
    protected:
      std::shared_ptr<BitArray> a;     /**< Bitarray the Reader is reading*/
      WType w;
      unsigned ix;        /**< Current index of Reader */ 
    public:
      Reader(const std::shared_ptr<BitArray> host, unsigned ix0=0);   /**< Constructor */
      unsigned index() const { return ix; }              /**< Return current index */
      template<unsigned n>
      WType read();                                          /**< Read n values */
      void seek(unsigned ix);                            /**< Search for ix in array */
    };
    
    /** Class for writing the bit array.
     *  Uses seek method from Reader.
     *  Note: writing off the end is undefined! */
    class Writer: public Reader {
    private:
      using Reader::seek;
    public:
      Writer(std::shared_ptr<BitArray> host, unsigned ix0=0);       /**< Default constructor */
      ~Writer();                                       
      template<unsigned n>
      void write(WType x);                                 /**< Write to array */
      void flush();                                    /**< Flush buffer */
    };
  };

  class FastBitArray {
    static const unsigned logc = 9, bitc = 1u<<logc;
    std::shared_ptr<BitArray > bits;
    std::vector<unsigned> chks;
  public:
    FastBitArray(unsigned len);
    class Builder {
      std::shared_ptr<FastBitArray > ref;
      typename BitArray::Writer w;
      unsigned chkpop;
      std::vector<unsigned> pchk;
    public:
      Builder(unsigned len);
      unsigned index() const { return w.index(); }
      template<unsigned n>
      void write(WType x);
      std::shared_ptr<FastBitArray > finish();
    };
  public:
    std::shared_ptr<BitArray> bit_array() const { return bits; }
    unsigned length() const { return bits->length(); }
    bool get(unsigned ix) const { return bits->get(ix); }
    unsigned count(unsigned ix0, unsigned ix1) const;
    unsigned find(unsigned ix0, unsigned nth) const;
  };
}
#endif
