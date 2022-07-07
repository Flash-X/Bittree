#include "Bittree_BitArray.h"
#include "Bittree_Bits.h"

namespace bittree {

  /** Constructor. Makes one extra word */
  BitArray::BitArray(unsigned len)
    : len_(len),
      wbuf_( (len+bitw)>>logw ) {
  }

  /**< get value of bit ix */
  bool BitArray::get(unsigned ix) const {
    // wbuf_ is our array of words.
    // Some tricks:
    //   a>>b == a/pow(2,b) and
    //   a<<b == a*pow(2,b)
    // So remember, when you see "1<<x", think "pow(2,x)".
    //
    // Therefor "ix>>logw" is a trick for quickly computing "ix/pow(2,logw)", and
    // is equivalent to "ix/bitw", which just gives us the index of the
    // word holding bit "ix".
    //
    // "ix & (bitw-1)" uses a different trick:
    //   a & (b-1) == a % b (works only when b is a power of 2)
    //
    // So, we're getting the word holding bit "ix" (ix/bitw), then shifting it down
    // by its position ordinal in the word (ix%bitw), then AND'ing by 1 to isolate just
    // that bit (zero's higher ones).
    return ix < len_ ? 1 & (wbuf_[ix>>logw] >> (ix & (bitw-1))) : 0;
  }

  /**< set value of bit ix */
  bool BitArray::set(unsigned ix, bool x) {
    // w0 means old word value, w1 is new word value
    WType w0 = wbuf_[ix>>logw];
    WType z = x ? WType(0) : ones;
    // We have w0, and we want to set bit at "ix%bitw" to have value "x".
    // Let's pretend x=true, then we want to set the bit to 1, which can be
    // accomiplished with:
    //   w1 = w0 | 1<<(ix%bitw).
    // Since x=true, z must be zero, so the actual code below reduces to the
    // previous line thanks to XOR'ing with 0 being a no-op.
    //
    // Now the case where x=false, then z=<all-ones>. XOR'ing against all ones
    // is just a bitflilp. So what the epxression below does is flip the bits into
    // "assume x=true" case, set the bit in the previous way (OR with a 1), and
    // then flip the bits back. Thus we have accomplished a set to 0 instead of
    // a set to 1.
    //
    // The utilty of this compact expression is that it handles both cases
    // without inducing if-else branch instructions.
    WType w1 = z^((z^w0) | one<<(ix&(bitw-1)));
    wbuf_[ix>>logw] = w1;
    return w1 != w0;
  }

  /** count 1's in whole array */
  unsigned BitArray::count() const {
    return count(0, len_);
  }

  /** count 1's in interval [ix0,ix1)
   * \todo add protection in case of out of bounds */
  unsigned BitArray::count(unsigned ix0, unsigned ix1) const {
    if(ix1 <= ix0) return 0;
    unsigned iw0 = ix0 >> logw;
    unsigned iw1 = (std::min(ix1, len_)-1) >> logw;
    WType m = ones << (ix0 & (bitw-1));
    unsigned pop = 0;
    for(unsigned iw=iw0; iw <= iw1; iw++) {
      if(iw == iw1)
        m &= ones >> (bitw-1-((std::min(ix1,len_)-1)&(bitw-1)));
      pop += static_cast<unsigned>(bitpop(wbuf_[iw] & m));
      m = ones;
    }
    return pop;
  }

  /** count 1's in either a or b */
  unsigned BitArray::count_xor(const BitArray& a, const BitArray& b,
                               unsigned ix0, unsigned ix1) {
    if(ix1 <= ix0) return 0;
    unsigned a_ix1 = std::min(ix1, a.len_);
    unsigned b_ix1 = std::min(ix1, b.len_);
    unsigned a_iw1 = (a_ix1 + bitw-1) >> logw;
    unsigned b_iw1 = (b_ix1 + bitw-1) >> logw;
    unsigned iw0 = ix0 >> logw;
    unsigned iw1 = std::max(a_iw1, b_iw1);
    WType m = ones << (ix0 & (bitw-1));
    unsigned pop = 0;
    for(unsigned iw=iw0; iw < iw1; iw++) {
      WType aw = iw < a_iw1 ? a.wbuf_[iw] : WType(0);
      aw &= ones >> (iw+1 < a_iw1 ? 0 : bitw-1-((a_ix1-1)&(bitw-1)));
      WType bw = iw < b_iw1 ? b.wbuf_[iw] : WType(0);
      bw &= ones >> (iw+1 < b_iw1 ? 0 : bitw-1-((b_ix1-1)&(bitw-1)));
      pop += static_cast<unsigned>(bitpop(m & (aw ^ bw)));
      m = ones;
    }
    return pop;
  }

  unsigned BitArray::find(unsigned ix0, unsigned nth) const {
    unsigned iw = ix0 >> logw;
    WType m = ones << (ix0&(bitw-1));
    while(true) {
      if(iw >= len_>>logw)
        m &= ones >> (bitw-1-((len_-1)&(bitw-1)));
      WType w = m & wbuf_[iw];
      unsigned pop = static_cast<unsigned>(bitpop(w));
      if(pop > nth) {
        while(nth--)
          w &= w - 1;
        return (iw<<logw) + static_cast<unsigned>(bitffs(w)) -1u;
      }
      nth -= pop;
      m = ones;
      iw += 1;
    }
  }

  /** Fill whole Bit Array */
  void BitArray::fill(bool x) {
    fill(x, 0, len_);
  }

  /** Fill part of Bit Array
   * \todo bounds check*/
  void BitArray::fill(bool x, unsigned ix0, unsigned ix1) {
    if(ix0 >= ix1) return;
    unsigned iw0 = ix0 >> logw, iw1 = (ix1-1) >> logw;
    WType z = x ? WType(0) : ones;
    WType m = ones << (ix0 & (bitw-1));
    for(unsigned iw=iw0; iw <= iw1; iw++) {
      if(iw == iw1)
        m &= ones >> (bitw-1-((ix1-1)&(bitw-1)));
      wbuf_[iw] = z ^ ((z ^ wbuf_[iw]) | m);
      m = ones;
    }
  }

  /** Constructor */
  BitArray::Reader::Reader(std::shared_ptr<const BitArray> host, unsigned ix0):
    a_(host),
    w_(host->wbuf_[0]),
    ix_(0) {
    seek(ix0);
  }

  /** Read n values */
  template<unsigned n>
  BitArray::WType BitArray::Reader::read() {
    WType ans;
    if(((ix_+n)&(bitw-1u)) > (ix_&(bitw-1u)))
      ans = (w_>>(ix_&(bitw-1u))) & ((one<<n)-1u);
    else {
      ans = w_>>(ix_&(bitw-1u));
      w_ = ((ix_+n)&~(bitw-1u)) < a_->len_ ? a_->wbuf_[(ix_>>logw)+1] : WType(0);
      ans |= (w_ & ((one<<((ix_+n)&(bitw-1u)))-1u)) << (bitw-(ix_&(bitw-1u)));
    }
    ix_ += n;
    return ans;
  }
 
  /** Search for ix in array */
  void BitArray::Reader::seek(unsigned ix) {
    w_ = a_->wbuf_[ix>>logw] ;
    ix_ = ix;
  }

  /** Constructor */
  BitArray::Writer::Writer(std::shared_ptr<BitArray> host, unsigned ix0)
   :a_(host),
    w_(host->wbuf_[0]),
    ix_(0) {
    seek(ix0);
  }

  /** Search for ix in array */
  void BitArray::Writer::seek(unsigned ix) {
    w_ = a_->wbuf_[ix>>logw] ;
    ix_ = ix;
  }

  /** Destructor */
  BitArray::Writer::~Writer() {
    flush();
  }

  /** Write n values to array
   * \todo Assert ix_ + n <= a_->length() */
  template<unsigned n>
  void BitArray::Writer::write(BitArray::WType x) {
    if(((ix_+n)&(bitw-1u)) > (ix_&(bitw-1u))) {
      WType m = ((one<<n)-one)<<(ix_&(bitw-1u));
      w_ = (w_ & ~m) | x<<(ix_&(bitw-1u));
    }
    else {
      WType m = ones<<(ix_&(bitw-1u));
      w_ = (w_ & ~m) | x<<(ix_&(bitw-1u));
      a_->wbuf_[ix_>>logw] = w_;
      w_ = a_->wbuf_[(ix_>>logw)+1];
      m = (one<<((ix_+n)&(bitw-1u)))-one;
      w_ = (w_ & ~m) | x>>(bitw-(ix_&(bitw-1u)));
    }
    ix_ += n;
  }

  /** Flush buffer */
  void BitArray::Writer::flush() {
    a_->wbuf_[ix_>>logw] = w_;
  }

  /** Constructor for FastBitArray.
    */
  FastBitArray::FastBitArray(unsigned len):
    BitArray(len),
    chks_(len>>logc) {
  }

  FastBitArray::Builder::Builder(unsigned len):
    a_(std::make_shared<FastBitArray>(len)),
    w_(BitArray::Writer(a_, 0)),
    chkpop_(0),
    pchk_(a_->chks_.data()) {
  }

  /** Wraps BitArray::Writer::write<n>, but also tracks cumulative
    * bitpop so FastBitArray can store it at intervals of bitc (=2^9).
    *
    * \todo combine branching statements if((ix&(bitc-1u))+n >= bitc)
    * \todo isolate n bits? (in case x is >= 2^n) )
    */
  template<unsigned n>
  void FastBitArray::Builder::write(WType x) {
    unsigned ix = w_.index();
    // NOTE: ix&(bitc-1u) = ix mod bitc
    // So this checks if the current index is within n bits
    //   of a checkpoint (where ix mod bitc = 0, ix/bitc = k).
    if((ix&(bitc-1u))+n >= bitc) {
      // Store cumulative bitpop up to bitc*2^k
      if(n == 1)
        *pchk_++ = ( chkpop_ + static_cast<unsigned>(x) );
      else {
        WType m = ~WType(0) >> (bitw-(bitc-(ix&(bitc-1u))));
        *pchk_++ = ( chkpop_ + static_cast<unsigned>(bitpop(x & m)) );
      }
    }
    // Keep track of total bitpop so far
    chkpop_ += n == 1 ? x : static_cast<unsigned>(bitpop(x));
    // Actually write bits
    w_.write<n>(x);
  }

  std::shared_ptr<FastBitArray> FastBitArray::Builder::finish() {
    w_.flush();
    return a_;
  }

  /** A theoretically faster algorithm for count making use of the cached cumulative pops
    * \todo Performance testing
    */
  unsigned FastBitArray::count(unsigned ix0, unsigned ix1) const {
    // If ix0 and ix1 are separated by over bitc bits, can use the chks_ array
    // to shorten the length to count.
    if(ix1>>logc > ix0>>logc) {
      // pop0 = total count up to next multiple of bitc (inclusive)
      // pop1 = total count up to previous multiple of bitc (exclusive)
      unsigned pop0 = ix0 == 0 ? 0 : chks_[((ix0+bitc-1u)>>logc)-1];
      unsigned pop1 = chks_[(ix1>>logc)-1];
      // (ix0+bitc-1u) & ~(bitc-1u) = the next multiple of bitc (inclusive)
      // ix1 & ~(bitc-1u) = the prevoius multiple of bitc (exclusive)
      return
        BitArray::count(ix0, (ix0+bitc-1u) & ~(bitc-1u)) +
        (pop1 - pop0) +
        BitArray::count(ix1 & ~(bitc-1u), ix1);
    }
    else
      return BitArray::count(ix0, ix1);
  }

  unsigned FastBitArray::find(unsigned ix0, unsigned nth) const {
    unsigned pop0 = count(0, ix0);
    unsigned a = (ix0+bitc-1u)>>logc;
    unsigned b = len_>>logc;
    if(a <= b && (a==0 ? 0 : chks_[a-1]) - pop0 <= nth) {
      // binary search of interval [a,b] (both inclusive)
      while(a + 10 < b) {
        unsigned x = (a+b)>>1;
        if(pop0 + nth < chks_[x-1]) // x cant be 0
          b = x-1;
        else
          a = x;
      }
      // linear search
      while(a < b && chks_[a]-pop0 <= nth)
        a += 1;
      nth -= (a==0 ? 0 : chks_[a-1]) - pop0;
      return BitArray::find(a<<logc, nth);
    }
    else
      return BitArray::find(ix0, nth);
  }

  template void FastBitArray::Builder::write<1u>(BitArray::WType x);
  template void FastBitArray::Builder::write<2u>(BitArray::WType x);
  template void FastBitArray::Builder::write<4u>(BitArray::WType x);
  template void FastBitArray::Builder::write<8u>(BitArray::WType x);
  template BitArray::WType BitArray::Reader::read<1u>();
  template BitArray::WType BitArray::Reader::read<2u>();
  template BitArray::WType BitArray::Reader::read<4u>();
  template BitArray::WType BitArray::Reader::read<8u>();

}
