#ifndef ANALYSISDATAFORMATS_BOOSTEDOBJECTS_MET_HH
#define ANALYSISDATAFORMATS_BOOSTEDOBJECTS_MET_HH

#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

namespace vlq{
  class Met ;
  typedef std::vector<vlq::Met> MetCollection ;
  class Met : public vlq::Candidate {
    public:
      Met() : 
        FullPhi_       (-1000),  
        FullPt_        (-1000),  
        FullPx_        (-1000),  
        FullPy_        (-1000),  
        FulluncorPhi_  (-1000),  
        FulluncorPt_   (-1000),  
        FulluncorSumEt_(-1000),  
        NoHFPhi_       (-1000),  
        NoHFPt_        (-1000),  
        NoHFPx_        (-1000),  
        NoHFPy_        (-1000),  
        NoHFuncorPhi_  (-1000),  
        NoHFuncorPt_   (-1000),  
        NoHFuncorSumEt_(-1000)    
    {}
      Met (const vlq::Met& met) :
        vlq::Candidate(met),
        FullPhi_        (met.getFullPhi        () ),
        FullPt_         (met.getFullPt         () ),
        FullPx_         (met.getFullPx         () ),
        FullPy_         (met.getFullPy         () ),
        FulluncorPhi_   (met.getFulluncorPhi   () ),
        FulluncorPt_    (met.getFulluncorPt    () ),
        FulluncorSumEt_ (met.getFulluncorSumEt () ),
        NoHFPhi_        (met.getNoHFPhi        () ),
        NoHFPt_         (met.getNoHFPt         () ),
        NoHFPx_         (met.getNoHFPx         () ),
        NoHFPy_         (met.getNoHFPy         () ),
        NoHFuncorPhi_   (met.getNoHFuncorPhi   () ),
        NoHFuncorPt_    (met.getNoHFuncorPt    () ),
        NoHFuncorSumEt_ (met.getNoHFuncorSumEt () )
    { }
      ~Met () {}

      float getFullPhi        () const {return FullPhi_        ;}
      float getFullPt         () const {return FullPt_         ;}
      float getFullPx         () const {return FullPx_         ;}
      float getFullPy         () const {return FullPy_         ;}
      float getFulluncorPhi   () const {return FulluncorPhi_   ;}
      float getFulluncorPt    () const {return FulluncorPt_    ;}
      float getFulluncorSumEt () const {return FulluncorSumEt_ ;}
      float getNoHFPhi        () const {return NoHFPhi_        ;}
      float getNoHFPt         () const {return NoHFPt_         ;}
      float getNoHFPx         () const {return NoHFPx_         ;}
      float getNoHFPy         () const {return NoHFPy_         ;}
      float getNoHFuncorPhi   () const {return NoHFuncorPhi_   ;}
      float getNoHFuncorPt    () const {return NoHFuncorPt_    ;}
      float getNoHFuncorSumEt () const {return NoHFuncorSumEt_ ;}

      void setFullPhi        ( const float& FullPhi )       { FullPhi_ = FullPhi               ;}
      void setFullPt         ( const float& FullPt )        { FullPt_ = FullPt                 ;}
      void setFullPx         ( const float& FullPx )        { FullPx_ = FullPx                 ;}
      void setFullPy         ( const float& FullPy )        { FullPy_ = FullPy                 ;}
      void setFulluncorPhi   ( const float& FulluncorPhi )  { FulluncorPhi_ = FulluncorPhi     ;}
      void setFulluncorPt    ( const float& FulluncorPt )   { FulluncorPt_ = FulluncorPt       ;}
      void setFulluncorSumEt ( const float& FulluncorSumEt ){ FulluncorSumEt_ = FulluncorSumEt ;}
      void setNoHFPhi        ( const float& NoHFPhi )       { NoHFPhi_ = NoHFPhi               ;}
      void setNoHFPt         ( const float& NoHFPt )        { NoHFPt_ = NoHFPt                 ;}
      void setNoHFPx         ( const float& NoHFPx )        { NoHFPx_ = NoHFPx                 ;}
      void setNoHFPy         ( const float& NoHFPy )        { NoHFPy_ = NoHFPy                 ;}
      void setNoHFuncorPhi   ( const float& NoHFuncorPhi )  { NoHFuncorPhi_ = NoHFuncorPhi     ;}
      void setNoHFuncorPt    ( const float& NoHFuncorPt )   { NoHFuncorPt_ = NoHFuncorPt       ;}
      void setNoHFuncorSumEt ( const float& NoHFuncorSumEt ){ NoHFuncorSumEt_ = NoHFuncorSumEt ;}

    protected:
      float FullPhi_        ;
      float FullPt_         ;
      float FullPx_         ;
      float FullPy_         ;
      float FulluncorPhi_   ;
      float FulluncorPt_    ;
      float FulluncorSumEt_ ;
      float NoHFPhi_        ;
      float NoHFPt_         ;
      float NoHFPx_         ;
      float NoHFPy_         ;
      float NoHFuncorPhi_   ;
      float NoHFuncorPt_    ;
      float NoHFuncorSumEt_ ;
  };
}

#endif
