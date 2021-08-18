#ifndef BITTREE_BITARRAY_H__
#define BITTREE_BITARRAY_H__

#include "Bittree_prelude.h"

namespace bittree {



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
    /** Definition of Word Type */
    typedef unsigned int WType;

    // Static variables
    /** Log_2 of memory allocated to objects of class W, in bits */
    static const unsigned logw = Log<2,sizeof(WType)*CHAR_BIT>::val; 
    static const unsigned bitw = 1u << logw; /**< bitwidth of WType */
    static const WType one = WType(1);    /**< 1 cast as WType */
    static const WType ones = ~WType(0);  /**< Maximum length string of binary 1s cast as WType */

  public:
    // Constructor
    BitArray(unsigned len);

    // Getters and setters
    unsigned length() const { return len_; }
    unsigned word_count() const { return (len_+bitw-1u)>>logw; }
    std::vector<WType> word_buf() { return wbuf_; }
    
    bool get(unsigned ix) const;
    bool set(unsigned ix, bool x);
    
    unsigned count(unsigned ix0, unsigned ix1) const;
    unsigned count() const;
    static unsigned count_xor(const BitArray& a, const BitArray& b,
                              unsigned ix0, unsigned ix1);
    
    unsigned find(unsigned ix0, unsigned nth) const;
    
    void fill(bool x);
    void fill(bool x, unsigned ix0, unsigned ix1);

  private:
    // Private members
    unsigned           len_;  /**< Length of Bit Array. Access with length() */
    std::vector<WType> wbuf_; /**< Word buffer of type WType. Access with word_buf() */

  public:

    /** Class for reading off the bit array.
     *  Note: reading off the end is safe and returns 0 bits   */
    class Reader {
    public:
      Reader(std::shared_ptr<const BitArray> host, unsigned ix0=0);
      unsigned index() const { return ix_; }  /**< Return current index */
      template<unsigned n>
      WType read();
      void seek(unsigned ix);
    protected:
      std::shared_ptr<const BitArray> a_; /**< Bitarray the Reader is reading*/
      WType w_;                           /**< Current word of Reader */
      unsigned ix_;                       /**< Current index of Reader */
    };
    
    /** Class for writing the bit array.
     *  Note: writing off the end is undefined! */
    class Writer {
    public:
      Writer(std::shared_ptr<BitArray> host, unsigned ix0=0);
      ~Writer();
      unsigned index() const { return ix_; }  /**< Return current index */
      template<unsigned n>
      void write(WType x);
      void seek(unsigned ix);
      void flush();
    protected:
      std::shared_ptr<BitArray> a_; /**< Bitarray the Writer is writing on*/
      WType w_;                     /**< Current word of Writer */
      unsigned ix_;                 /**< Current index of Writer */
    };
  };

  /** FastBitArray class
    *
    */
  class FastBitArray {

    static const unsigned logc = 9;
    static const unsigned bitc = 1u<<logc;

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
      void write(BitArray::WType x);
      std::shared_ptr<FastBitArray > finish();
    };

    std::shared_ptr<BitArray> bit_array() const { return bits_; }
    unsigned length() const { return bits_->length(); }
    bool get(unsigned ix) const { return bits_->get(ix); }
    unsigned count(unsigned ix0, unsigned ix1) const;
    unsigned find(unsigned ix0, unsigned nth) const;

  private:
    std::shared_ptr<BitArray > bits_;
    std::vector<unsigned> chks_;
  };
}
#endif
