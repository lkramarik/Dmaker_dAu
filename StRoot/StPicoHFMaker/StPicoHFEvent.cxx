#include "StPicoEvent/StPicoEvent.h"
#include "StPicoHFEvent.h"
#include "StHFPair.h"

ClassImp(StPicoHFEvent)

TClonesArray *StPicoHFEvent::fgHFSecondaryVerticesArray = 0;
// _________________________________________________________
StPicoHFEvent::StPicoHFEvent() : mRunId(-1), mEventId(-1), mNHFSecondaryVertices(0),
						  mHFSecondaryVerticesArray(NULL) {
  // -- Default constructor
  if (!fgHFSecondaryVerticesArray) fgHFSecondaryVerticesArray = new TClonesArray("StHFPair");
  mHFSecondaryVerticesArray = fgHFSecondaryVerticesArray;
}

// _________________________________________________________
StPicoHFEvent::StPicoHFEvent(unsigned int mode) : mRunId(-1), mEventId(-1), mNHFSecondaryVertices(0),
						  mHFSecondaryVerticesArray(NULL) {
  // -- Constructor with mode selection
  if (mode == StPicoHFEvent::kTwoAndTwoParticleDecay) {
    if (!fgHFSecondaryVerticesArray) fgHFSecondaryVerticesArray = new TClonesArray("StHFPair");
    mHFSecondaryVerticesArray = fgHFSecondaryVerticesArray;
  }
  else if (mode == StPicoHFEvent::kTwoParticleDecay) {
    if (!fgHFSecondaryVerticesArray) fgHFSecondaryVerticesArray = new TClonesArray("StHFPair");
    mHFSecondaryVerticesArray = fgHFSecondaryVerticesArray;
  }
  else {
    if (!fgHFSecondaryVerticesArray) fgHFSecondaryVerticesArray = new TClonesArray("StHFPair");
    mHFSecondaryVerticesArray = fgHFSecondaryVerticesArray;
  }
}

// _________________________________________________________
void StPicoHFEvent::addPicoEvent(StPicoEvent const & picoEvent) {
   // -- add StPicoEvent variables
   mRunId   = picoEvent.runId();
   mEventId = picoEvent.eventId();
}

// _________________________________________________________
void StPicoHFEvent::clear(char const *option) {
  mHFSecondaryVerticesArray->Clear(option);
  mRunId                = -1;
  mEventId              = -1;
  mNHFSecondaryVertices = 0;
}

// _________________________________________________________
void StPicoHFEvent::addHFSecondaryVertexPair(StHFPair const* t) {
  TClonesArray &vertexArray = *mHFSecondaryVerticesArray;
  new(vertexArray[mNHFSecondaryVertices++]) StHFPair(t);
}