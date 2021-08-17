#include "Bittree_bitarray.h"
#include "Bittree_bits.h"
#include "Bittree_mem.h"
#include "Bittree_ref.h"

namespace BitTree {
  Ref<BitArray > BitArray::make(unsigned n) {
    Ref<BitArray > ref;
    unsigned nw = (n + bitw)>>logw; // one extra
    Mem mem = {
      alignof(BitArray),
      sizeof(BitArray) + (nw-1)*sizeof(WType)
    };
    BitArray *it = new(ref.alloc(mem)) BitArray(n);
    it->wbuf[nw-1] = WType(0);
    return ref;
  }

  bool BitArray::get(unsigned ix) const {
    // wbuf is our array of words.
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
    return ix < this->len ? 1 & (wbuf[ix>>logw] >> (ix & (bitw-1))) : 0;
  }

  bool BitArray::set(unsigned ix, bool x) {
    // w0 means old word value, w1 is new word value
    WType w0 = wbuf[ix>>logw];
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
    wbuf[ix>>logw] = w1;
    return w1 != w0;
  }

  unsigned BitArray::count() const {
    return this->count(0, this->len);
  }

  unsigned BitArray::count(unsigned ix0, unsigned ix1) const {
    if(ix1 <= ix0) return 0;
    unsigned iw0 = ix0 >> logw;
    unsigned iw1 = (std::min(ix1, this->len)-1) >> logw;
    WType m = ones << (ix0 & (bitw-1));
    unsigned pop = 0;
    for(unsigned iw=iw0; iw <= iw1; iw++) {
      if(iw == iw1)
        m &= ones >> (bitw-1-((std::min(ix1,this->len)-1)&(bitw-1)));
      pop += static_cast<unsigned>(bitpop(wbuf[iw] & m));
      m = ones;
    }
    return pop;
  }

  unsigned BitArray::count_xor(
      const BitArray *a,
      const BitArray *b,
      unsigned ix0, unsigned ix1
    ) {
    if(ix1 <= ix0) return 0;
    unsigned a_ix1 = std::min(ix1, a->len);
    unsigned b_ix1 = std::min(ix1, b->len);
    unsigned a_iw1 = (a_ix1 + bitw-1) >> logw;
    unsigned b_iw1 = (b_ix1 + bitw-1) >> logw;
    unsigned iw0 = ix0 >> logw;
    unsigned iw1 = std::max(a_iw1, b_iw1);
    WType m = ones << (ix0 & (bitw-1));
    unsigned pop = 0;
    for(unsigned iw=iw0; iw < iw1; iw++) {
      WType aw = iw < a_iw1 ? a->wbuf[iw] : WType(0);
      aw &= ones >> (iw+1 < a_iw1 ? 0 : bitw-1-((a_ix1-1)&(bitw-1)));
      WType bw = iw < b_iw1 ? b->wbuf[iw] : WType(0);
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
      if(iw >= len>>logw)
        m &= ones >> (bitw-1-((len-1)&(bitw-1)));
      WType w = m & wbuf[iw];
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

  void BitArray::fill(bool x) {
    this->fill(x, 0, this->len);
  }

  void BitArray::fill(bool x, unsigned ix0, unsigned ix1) {
    if(ix0 >= ix1) return;
    DBG_ASSERT(ix1 <= this->len);
    unsigned iw0 = ix0 >> logw, iw1 = (ix1-1) >> logw;
    WType z = x ? WType(0) : ones;
    WType m = ones << (ix0 & (bitw-1));
    for(unsigned iw=iw0; iw <= iw1; iw++) {
      if(iw == iw1)
        m &= ones >> (bitw-1-((ix1-1)&(bitw-1)));
      wbuf[iw] = z ^ ((z ^ wbuf[iw]) | m);
      m = ones;
    }
  }

  BitArray::Reader::Reader(const BitArray *a, unsigned ix0):
    a(const_cast<BitArray*>(a)),
    w( a->wbuf[0] ),
    ix(0) {
    this->seek(ix0);
  }

  template<unsigned n>
  WType BitArray::Reader::read() {
    unsigned n_ = n;
    BitArray *a = this->a;
    WType &w = this->w;
    unsigned &ix = this->ix;
    WType ans;
    if(((ix+n)&(bitw-1u)) > (ix&(bitw-1u)))
      ans = (w>>(ix&(bitw-1u))) & ((one<<n_)-1u);
    else {
      ans = w>>(ix&(bitw-1u));
      w = ((ix+n)&~(bitw-1u)) < a->len ? a->wbuf[(ix>>logw)+1] : WType(0);
      ans |= (w & ((one<<((ix+n)&(bitw-1u)))-1u)) << (bitw-(ix&(bitw-1u)));
    }
    ix += n;
    return ans;
  }
  
  void BitArray::Reader::seek(unsigned ix) {
    this->w = a->wbuf[ix>>logw] ;
    this->ix = ix;
  }

  BitArray::Writer::Writer(BitArray *host, unsigned ix0):
    BitArray::Reader(host, ix0) {
  }
  BitArray::Writer::~Writer() {
    this->flush();
  }

  template<unsigned n>
  void BitArray::Writer::write(WType x) {
    unsigned n_ = n;
    BitArray *a = this->a;
    WType &w = this->w;
    unsigned &ix = this->ix;
    DBG_ASSERT(ix + n <= a->length());
    if(((ix+n)&(bitw-1u)) > (ix&(bitw-1u))) {
      WType m = ((one<<n_)-one)<<(ix&(bitw-1u));
      w = (w & ~m) | x<<(ix&(bitw-1u));
    }
    else {
      WType m = ones<<(ix&(bitw-1u));
      w = (w & ~m) | x<<(ix&(bitw-1u));
      a->wbuf[ix>>logw] = w;
      w = a->wbuf[(ix>>logw)+1];
      m = (one<<((ix+n)&(bitw-1u)))-one;
      w = (w & ~m) | x>>(bitw-(ix&(bitw-1u)));
    }
    ix += n;
  }

  void BitArray::Writer::flush() {
    this->a->wbuf[this->ix>>logw] = this->w;
  }

  FastBitArray::FastBitArray(unsigned len):
    bitsref(BitArray::make(len)),
    bits(bitsref),
    chksref(Ref_::new_array<unsigned>(len>>logc)),
    chks(chksref) {
  }

  FastBitArray::Builder::Builder(unsigned len):
    ref(Ref_::new1<FastBitArray,unsigned>(len)),
    w(typename BitArray::Writer(ref->bits, 0)),
    pchk(ref->chks),
    chkpop(0) {
  }

  template<unsigned n>
  void FastBitArray::Builder::write(WType x) {
    unsigned ix = w.index();
    if((ix&(bitc-1u))+n >= bitc) {
      if(n == 1)
        *pchk++ = chkpop + static_cast<unsigned>(x);
      else {
        WType m = ~WType(0) >> (BitArray::bitw-(bitc-(ix&(bitc-1u))));
        *pchk++ = chkpop + static_cast<unsigned>(bitpop(x & m));
      }
    }
    chkpop += n == 1 ? x : static_cast<unsigned>(bitpop(x));
    w.write<n>(x);
  }

  Ref<FastBitArray > FastBitArray::Builder::finish() {
    w.flush();
    return ref;
  }

  unsigned FastBitArray::count(unsigned ix0, unsigned ix1) const {
    if(ix1>>logc > ix0>>logc) {
      unsigned pop0 = ix0 == 0 ? 0 : chks[((ix0+bitc-1u)>>logc)-1];
      unsigned pop1 = chks[(ix1>>logc)-1];
      return
        bits->count(ix0, (ix0+bitc-1u) & ~(bitc-1u)) +
        (pop1 - pop0) +
        bits->count(ix1 & ~(bitc-1u), ix1);
    }
    else
      return bits->count(ix0, ix1);
  }

  unsigned FastBitArray::find(unsigned ix0, unsigned nth) const {
    unsigned pop0 = this->count(0, ix0);
    unsigned a = (ix0+bitc-1u)>>logc;
    unsigned b = bits->length()>>logc;
    if(a <= b && (a==0 ? 0 : chks[a-1]) - pop0 <= nth) {
      // binary search of interval [a,b] (both inclusive)
      while(a + 10 < b) {
        unsigned x = (a+b)>>1;
        if(pop0 + nth < chks[x-1]) // x cant be 0
          b = x-1;
        else
          a = x;
      }
      // linear search
      while(a < b && chks[a]-pop0 <= nth)
        a += 1;
      nth -= (a==0 ? 0 : chks[a-1]) - pop0;
      return bits->find(a<<logc, nth);
    }
    else
      return bits->find(ix0, nth);
  }

  template void FastBitArray::Builder::write<1u>(WType x);
  template void FastBitArray::Builder::write<2u>(WType x);
  template void FastBitArray::Builder::write<4u>(WType x);
  template void FastBitArray::Builder::write<8u>(WType x);
  template WType BitArray::Reader::read<1u>();
  template WType BitArray::Reader::read<2u>();
  template WType BitArray::Reader::read<4u>();
  template WType BitArray::Reader::read<8u>();

}
