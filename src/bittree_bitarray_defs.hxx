#ifndef _1345010d_9317_43fa_90c7_d2dc7bc16926
#define _1345010d_9317_43fa_90c7_d2dc7bc16926

#include "bittree_prelude.hxx"
#include "bittree_ref_defs.hxx"

namespace BitTree {
  /** Stores, reads, and writes Bit Arrays.
   *    */
  template<class W>
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
    static const unsigned logw = Log<2,sizeof(W)*CHAR_BIT>::val; 
    /** bitwidth of W */
    static const unsigned bitw = 1u << logw;
    static const W one = W(1);    /**< 1 cast as W */
    static const W ones = ~W(0);  /**< Maximum length string of binary 1s cast as W */
  private:
    unsigned len;       /**< Length of Bit Array. Access with length() */
    W wbuf[1];          /**< Word buffer of type W. Access with word_buf() */
  private:
    BitArray(unsigned len): len(len) {}
  public:
    static Ref<BitArray<W> > make(unsigned n);
  public:
    unsigned length() const { return len; }
    unsigned word_count() const { return (len+bitw-1u)>>logw; }
    W* word_buf() { return wbuf; }
    
    bool get(unsigned ix) const;     /**< get value of bit ix */
    bool set(unsigned ix, bool x);   /**< set value of bit ix */
    
    /** count 1's in interval [ix0,ix1) */
    unsigned count(unsigned ix0, unsigned ix1) const;
    /** count 1's in whole array */
    unsigned count() const;
    /** count 1's in either a or b */
    static unsigned count_xor(
      const BitArray<W> *a,
      const BitArray<W> *b,
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
      BitArray<W> *a;     /**< Bitarray the Reader is reading*/
      W w;
      unsigned ix;        /**< Current index of Reader */ 
    public:
      Reader(const BitArray<W> *host, unsigned ix0=0);   /**< Constructor */
      unsigned index() const { return ix; }              /**< Return current index */
      template<unsigned n>
      W read();                                          /**< Read n values */
      void seek(unsigned ix);                            /**< Search for ix in array */
    };
    
    /** Class for writing the bit array.
     *  Uses seek method from Reader.
     *  Note: writing off the end is undefined! */
    class Writer: public Reader {
    private:
      using Reader::seek;
    public:
      Writer(BitArray<W> *host, unsigned ix0=0);       /**< Default constructor */
      ~Writer();                                       
      template<unsigned n>
      void write(W x);                                 /**< Write to array */
      void flush();                                    /**< Flush buffer */
    };
  };

  template<class W>
  class FastBitArray {
    static const unsigned logc = 9, bitc = 1u<<logc;
    Ref<BitArray<W> > bitsref;
    BitArray<W> *bits;
    Ref<unsigned> chksref;
    unsigned *chks;
  public:
    FastBitArray(unsigned len);
    class Builder {
      Ref<FastBitArray<W> > ref;
      typename BitArray<W>::Writer w;
      unsigned *pchk, chkpop;
    public:
      Builder(unsigned len);
      unsigned index() const { return w.index(); }
      template<unsigned n>
      void write(W x);
      Ref<FastBitArray<W> > finish();
    };
  public:
    const BitArray<W>* bit_array() const { return bits; }
    unsigned length() const { return bits->length(); }
    bool get(unsigned ix) const { return bits->get(ix); }
    unsigned count(unsigned ix0, unsigned ix1) const;
    unsigned find(unsigned ix0, unsigned nth) const;
  };
}
#endif
