///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TriVDCHit                                                                 //
//                                                                           //
// Class representing a single hit for the VDC                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TriVDCHit.h"
#include "TriVDCTimeToDistConv.h"

ClassImp(TriVDCHit)

const Double_t TriVDCHit::kBig = 1.e38; // Arbitrary large value

//_____________________________________________________________________________
Double_t TriVDCHit::ConvertTimeToDist(Double_t slope)
{
  // Converts TDC time to drift distance
  // Takes the (estimated) slope of the track as an argument
  
  TriVDCTimeToDistConv* ttdConv = (fWire) ? fWire->GetTTDConv() : NULL;
  
  if (ttdConv) {
    // If a time to distance algorithm exists, use it to convert the TDC time 
    // to the drift distance
    fDist = ttdConv->ConvertTimeToDist(fTime, slope, &fdDist);
    return fDist;
  }
  
  Error("ConvertTimeToDist()", "No Time to dist algorithm available");
  return 0.0;

}

//_____________________________________________________________________________
Int_t TriVDCHit::Compare( const TObject* obj ) const 
{
  // Used to sort hits
  // A hit is "less than" another hit if it occurred on a lower wire number.
  // Also, for hits on the same wire, the first hit on the wire (the one with
  // the smallest time) is "less than" one with a higher time.  If the hits
  // are sorted according to this scheme, they will be in order of increasing
  // wire number and, for each wire, will be in the order in which they hit
  // the wire

  if( !obj || IsA() != obj->IsA() || !fWire )
    return -1;

  const TriVDCHit* hit = static_cast<const TriVDCHit*>( obj );
 
  Int_t myWireNum = fWire->GetNum();
  Int_t hitWireNum = hit->GetWire()->GetNum();
  // Compare wire numbers
  if (myWireNum < hitWireNum) return -1;
  if (myWireNum > hitWireNum) return  1;
  if (myWireNum == hitWireNum) {
    // If wire numbers are the same, compare times
    Double_t hitTime = hit->GetTime();
    if (fTime < hitTime) return -1;
    if (fTime > hitTime) return  1;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
