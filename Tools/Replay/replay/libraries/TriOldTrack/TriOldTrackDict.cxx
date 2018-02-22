// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TriOldTrackDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TriHRS.h"
#include "TriSpectrometer.h"
#include "TriSpectrometerDetector.h"
#include "TriTrackingDetector.h"
#include "TriVDC.h"
#include "TriVDCUVPlane.h"
#include "TriVDCPlane.h"
#include "TriVDCCluster.h"
#include "TriVDCHit.h"
#include "TriVDCTimeToDistConv.h"
#include "TriVDCUVTrack.h"
#include "TriVDCTrackPair.h"
#include "TriVDCTrackID.h"
#include "TriVDCWire.h"
#include "TriTrack.h"
#include "TriTrackInfo.h"
#include "TriTrackingModule.h"
#include "TriNonTrackingDetector.h"
#include "TriScintillator.h"
#include "TriTrackID.h"
#include "TriVDCAnalyticTTDConv.h"
#include "TriTriggerTime.h"
#include "TriTrackProj.h"
#include "TriXscin.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TriTrackInfo(void *p = 0);
   static void *newArray_TriTrackInfo(Long_t size, void *p);
   static void delete_TriTrackInfo(void *p);
   static void deleteArray_TriTrackInfo(void *p);
   static void destruct_TriTrackInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriTrackInfo*)
   {
      ::TriTrackInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriTrackInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriTrackInfo", ::TriTrackInfo::Class_Version(), "TriTrackInfo.h", 14,
                  typeid(::TriTrackInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriTrackInfo::Dictionary, isa_proxy, 4,
                  sizeof(::TriTrackInfo) );
      instance.SetNew(&new_TriTrackInfo);
      instance.SetNewArray(&newArray_TriTrackInfo);
      instance.SetDelete(&delete_TriTrackInfo);
      instance.SetDeleteArray(&deleteArray_TriTrackInfo);
      instance.SetDestructor(&destruct_TriTrackInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriTrackInfo*)
   {
      return GenerateInitInstanceLocal((::TriTrackInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriTrackInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriTrackingModule(void *p = 0);
   static void *newArray_TriTrackingModule(Long_t size, void *p);
   static void delete_TriTrackingModule(void *p);
   static void deleteArray_TriTrackingModule(void *p);
   static void destruct_TriTrackingModule(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriTrackingModule*)
   {
      ::TriTrackingModule *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriTrackingModule >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriTrackingModule", ::TriTrackingModule::Class_Version(), "TriTrackingModule.h", 15,
                  typeid(::TriTrackingModule), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriTrackingModule::Dictionary, isa_proxy, 4,
                  sizeof(::TriTrackingModule) );
      instance.SetNew(&new_TriTrackingModule);
      instance.SetNewArray(&newArray_TriTrackingModule);
      instance.SetDelete(&delete_TriTrackingModule);
      instance.SetDeleteArray(&deleteArray_TriTrackingModule);
      instance.SetDestructor(&destruct_TriTrackingModule);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriTrackingModule*)
   {
      return GenerateInitInstanceLocal((::TriTrackingModule*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriTrackingModule*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriSpectrometer(void *p);
   static void deleteArray_TriSpectrometer(void *p);
   static void destruct_TriSpectrometer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriSpectrometer*)
   {
      ::TriSpectrometer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriSpectrometer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriSpectrometer", ::TriSpectrometer::Class_Version(), "TriSpectrometer.h", 23,
                  typeid(::TriSpectrometer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriSpectrometer::Dictionary, isa_proxy, 4,
                  sizeof(::TriSpectrometer) );
      instance.SetDelete(&delete_TriSpectrometer);
      instance.SetDeleteArray(&deleteArray_TriSpectrometer);
      instance.SetDestructor(&destruct_TriSpectrometer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriSpectrometer*)
   {
      return GenerateInitInstanceLocal((::TriSpectrometer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriSpectrometer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriHRS(void *p);
   static void deleteArray_TriHRS(void *p);
   static void destruct_TriHRS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriHRS*)
   {
      ::TriHRS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriHRS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriHRS", ::TriHRS::Class_Version(), "TriHRS.h", 14,
                  typeid(::TriHRS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriHRS::Dictionary, isa_proxy, 4,
                  sizeof(::TriHRS) );
      instance.SetDelete(&delete_TriHRS);
      instance.SetDeleteArray(&deleteArray_TriHRS);
      instance.SetDestructor(&destruct_TriHRS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriHRS*)
   {
      return GenerateInitInstanceLocal((::TriHRS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriHRS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriSpectrometerDetector(void *p);
   static void deleteArray_TriSpectrometerDetector(void *p);
   static void destruct_TriSpectrometerDetector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriSpectrometerDetector*)
   {
      ::TriSpectrometerDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriSpectrometerDetector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriSpectrometerDetector", ::TriSpectrometerDetector::Class_Version(), "TriSpectrometerDetector.h", 19,
                  typeid(::TriSpectrometerDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriSpectrometerDetector::Dictionary, isa_proxy, 4,
                  sizeof(::TriSpectrometerDetector) );
      instance.SetDelete(&delete_TriSpectrometerDetector);
      instance.SetDeleteArray(&deleteArray_TriSpectrometerDetector);
      instance.SetDestructor(&destruct_TriSpectrometerDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriSpectrometerDetector*)
   {
      return GenerateInitInstanceLocal((::TriSpectrometerDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriSpectrometerDetector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriTrackingDetector(void *p);
   static void deleteArray_TriTrackingDetector(void *p);
   static void destruct_TriTrackingDetector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriTrackingDetector*)
   {
      ::TriTrackingDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriTrackingDetector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriTrackingDetector", ::TriTrackingDetector::Class_Version(), "TriTrackingDetector.h", 16,
                  typeid(::TriTrackingDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriTrackingDetector::Dictionary, isa_proxy, 4,
                  sizeof(::TriTrackingDetector) );
      instance.SetDelete(&delete_TriTrackingDetector);
      instance.SetDeleteArray(&deleteArray_TriTrackingDetector);
      instance.SetDestructor(&destruct_TriTrackingDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriTrackingDetector*)
   {
      return GenerateInitInstanceLocal((::TriTrackingDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriTrackingDetector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriVDC(void *p);
   static void deleteArray_TriVDC(void *p);
   static void destruct_TriVDC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDC*)
   {
      ::TriVDC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDC", ::TriVDC::Class_Version(), "TriVDC.h", 18,
                  typeid(::TriVDC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDC::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDC) );
      instance.SetDelete(&delete_TriVDC);
      instance.SetDeleteArray(&deleteArray_TriVDC);
      instance.SetDestructor(&destruct_TriVDC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDC*)
   {
      return GenerateInitInstanceLocal((::TriVDC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCPlane(void *p = 0);
   static void *newArray_TriVDCPlane(Long_t size, void *p);
   static void delete_TriVDCPlane(void *p);
   static void deleteArray_TriVDCPlane(void *p);
   static void destruct_TriVDCPlane(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCPlane*)
   {
      ::TriVDCPlane *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCPlane >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCPlane", ::TriVDCPlane::Class_Version(), "TriVDCPlane.h", 22,
                  typeid(::TriVDCPlane), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCPlane::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCPlane) );
      instance.SetNew(&new_TriVDCPlane);
      instance.SetNewArray(&newArray_TriVDCPlane);
      instance.SetDelete(&delete_TriVDCPlane);
      instance.SetDeleteArray(&deleteArray_TriVDCPlane);
      instance.SetDestructor(&destruct_TriVDCPlane);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCPlane*)
   {
      return GenerateInitInstanceLocal((::TriVDCPlane*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCPlane*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCUVPlane(void *p = 0);
   static void *newArray_TriVDCUVPlane(Long_t size, void *p);
   static void delete_TriVDCUVPlane(void *p);
   static void deleteArray_TriVDCUVPlane(void *p);
   static void destruct_TriVDCUVPlane(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCUVPlane*)
   {
      ::TriVDCUVPlane *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCUVPlane >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCUVPlane", ::TriVDCUVPlane::Class_Version(), "TriVDCUVPlane.h", 18,
                  typeid(::TriVDCUVPlane), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCUVPlane::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCUVPlane) );
      instance.SetNew(&new_TriVDCUVPlane);
      instance.SetNewArray(&newArray_TriVDCUVPlane);
      instance.SetDelete(&delete_TriVDCUVPlane);
      instance.SetDeleteArray(&deleteArray_TriVDCUVPlane);
      instance.SetDestructor(&destruct_TriVDCUVPlane);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCUVPlane*)
   {
      return GenerateInitInstanceLocal((::TriVDCUVPlane*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCUVPlane*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCCluster(void *p = 0);
   static void *newArray_TriVDCCluster(Long_t size, void *p);
   static void delete_TriVDCCluster(void *p);
   static void deleteArray_TriVDCCluster(void *p);
   static void destruct_TriVDCCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCCluster*)
   {
      ::TriVDCCluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCCluster >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCCluster", ::TriVDCCluster::Class_Version(), "TriVDCCluster.h", 17,
                  typeid(::TriVDCCluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCCluster::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCCluster) );
      instance.SetNew(&new_TriVDCCluster);
      instance.SetNewArray(&newArray_TriVDCCluster);
      instance.SetDelete(&delete_TriVDCCluster);
      instance.SetDeleteArray(&deleteArray_TriVDCCluster);
      instance.SetDestructor(&destruct_TriVDCCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCCluster*)
   {
      return GenerateInitInstanceLocal((::TriVDCCluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCCluster*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCWire(void *p = 0);
   static void *newArray_TriVDCWire(Long_t size, void *p);
   static void delete_TriVDCWire(void *p);
   static void deleteArray_TriVDCWire(void *p);
   static void destruct_TriVDCWire(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCWire*)
   {
      ::TriVDCWire *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCWire >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCWire", ::TriVDCWire::Class_Version(), "TriVDCWire.h", 13,
                  typeid(::TriVDCWire), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCWire::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCWire) );
      instance.SetNew(&new_TriVDCWire);
      instance.SetNewArray(&newArray_TriVDCWire);
      instance.SetDelete(&delete_TriVDCWire);
      instance.SetDeleteArray(&deleteArray_TriVDCWire);
      instance.SetDestructor(&destruct_TriVDCWire);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCWire*)
   {
      return GenerateInitInstanceLocal((::TriVDCWire*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCWire*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCHit(void *p = 0);
   static void *newArray_TriVDCHit(Long_t size, void *p);
   static void delete_TriVDCHit(void *p);
   static void deleteArray_TriVDCHit(void *p);
   static void destruct_TriVDCHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCHit*)
   {
      ::TriVDCHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCHit", ::TriVDCHit::Class_Version(), "TriVDCHit.h", 14,
                  typeid(::TriVDCHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCHit::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCHit) );
      instance.SetNew(&new_TriVDCHit);
      instance.SetNewArray(&newArray_TriVDCHit);
      instance.SetDelete(&delete_TriVDCHit);
      instance.SetDeleteArray(&deleteArray_TriVDCHit);
      instance.SetDestructor(&destruct_TriVDCHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCHit*)
   {
      return GenerateInitInstanceLocal((::TriVDCHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriVDCTimeToDistConv(void *p);
   static void deleteArray_TriVDCTimeToDistConv(void *p);
   static void destruct_TriVDCTimeToDistConv(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCTimeToDistConv*)
   {
      ::TriVDCTimeToDistConv *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCTimeToDistConv >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCTimeToDistConv", ::TriVDCTimeToDistConv::Class_Version(), "TriVDCTimeToDistConv.h", 14,
                  typeid(::TriVDCTimeToDistConv), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCTimeToDistConv::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCTimeToDistConv) );
      instance.SetDelete(&delete_TriVDCTimeToDistConv);
      instance.SetDeleteArray(&deleteArray_TriVDCTimeToDistConv);
      instance.SetDestructor(&destruct_TriVDCTimeToDistConv);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCTimeToDistConv*)
   {
      return GenerateInitInstanceLocal((::TriVDCTimeToDistConv*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCTimeToDistConv*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCUVTrack(void *p = 0);
   static void *newArray_TriVDCUVTrack(Long_t size, void *p);
   static void delete_TriVDCUVTrack(void *p);
   static void deleteArray_TriVDCUVTrack(void *p);
   static void destruct_TriVDCUVTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCUVTrack*)
   {
      ::TriVDCUVTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCUVTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCUVTrack", ::TriVDCUVTrack::Class_Version(), "TriVDCUVTrack.h", 17,
                  typeid(::TriVDCUVTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCUVTrack::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCUVTrack) );
      instance.SetNew(&new_TriVDCUVTrack);
      instance.SetNewArray(&newArray_TriVDCUVTrack);
      instance.SetDelete(&delete_TriVDCUVTrack);
      instance.SetDeleteArray(&deleteArray_TriVDCUVTrack);
      instance.SetDestructor(&destruct_TriVDCUVTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCUVTrack*)
   {
      return GenerateInitInstanceLocal((::TriVDCUVTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCUVTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCTrackPair(void *p = 0);
   static void *newArray_TriVDCTrackPair(Long_t size, void *p);
   static void delete_TriVDCTrackPair(void *p);
   static void deleteArray_TriVDCTrackPair(void *p);
   static void destruct_TriVDCTrackPair(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCTrackPair*)
   {
      ::TriVDCTrackPair *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCTrackPair >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCTrackPair", ::TriVDCTrackPair::Class_Version(), "TriVDCTrackPair.h", 15,
                  typeid(::TriVDCTrackPair), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCTrackPair::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCTrackPair) );
      instance.SetNew(&new_TriVDCTrackPair);
      instance.SetNewArray(&newArray_TriVDCTrackPair);
      instance.SetDelete(&delete_TriVDCTrackPair);
      instance.SetDeleteArray(&deleteArray_TriVDCTrackPair);
      instance.SetDestructor(&destruct_TriVDCTrackPair);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCTrackPair*)
   {
      return GenerateInitInstanceLocal((::TriVDCTrackPair*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCTrackPair*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriTrackID(void *p);
   static void deleteArray_TriTrackID(void *p);
   static void destruct_TriTrackID(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriTrackID*)
   {
      ::TriTrackID *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriTrackID >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriTrackID", ::TriTrackID::Class_Version(), "TriTrackID.h", 13,
                  typeid(::TriTrackID), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriTrackID::Dictionary, isa_proxy, 4,
                  sizeof(::TriTrackID) );
      instance.SetDelete(&delete_TriTrackID);
      instance.SetDeleteArray(&deleteArray_TriTrackID);
      instance.SetDestructor(&destruct_TriTrackID);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriTrackID*)
   {
      return GenerateInitInstanceLocal((::TriTrackID*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriTrackID*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCTrackID(void *p = 0);
   static void *newArray_TriVDCTrackID(Long_t size, void *p);
   static void delete_TriVDCTrackID(void *p);
   static void deleteArray_TriVDCTrackID(void *p);
   static void destruct_TriVDCTrackID(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCTrackID*)
   {
      ::TriVDCTrackID *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCTrackID >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCTrackID", ::TriVDCTrackID::Class_Version(), "TriVDCTrackID.h", 13,
                  typeid(::TriVDCTrackID), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCTrackID::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCTrackID) );
      instance.SetNew(&new_TriVDCTrackID);
      instance.SetNewArray(&newArray_TriVDCTrackID);
      instance.SetDelete(&delete_TriVDCTrackID);
      instance.SetDeleteArray(&deleteArray_TriVDCTrackID);
      instance.SetDestructor(&destruct_TriVDCTrackID);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCTrackID*)
   {
      return GenerateInitInstanceLocal((::TriVDCTrackID*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCTrackID*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriTrack(void *p = 0);
   static void *newArray_TriTrack(Long_t size, void *p);
   static void delete_TriTrack(void *p);
   static void deleteArray_TriTrack(void *p);
   static void destruct_TriTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriTrack*)
   {
      ::TriTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriTrack", ::TriTrack::Class_Version(), "TriTrack.h", 20,
                  typeid(::TriTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriTrack::Dictionary, isa_proxy, 4,
                  sizeof(::TriTrack) );
      instance.SetNew(&new_TriTrack);
      instance.SetNewArray(&newArray_TriTrack);
      instance.SetDelete(&delete_TriTrack);
      instance.SetDeleteArray(&deleteArray_TriTrack);
      instance.SetDestructor(&destruct_TriTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriTrack*)
   {
      return GenerateInitInstanceLocal((::TriTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TriNonTrackingDetector(void *p);
   static void deleteArray_TriNonTrackingDetector(void *p);
   static void destruct_TriNonTrackingDetector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriNonTrackingDetector*)
   {
      ::TriNonTrackingDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriNonTrackingDetector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriNonTrackingDetector", ::TriNonTrackingDetector::Class_Version(), "TriNonTrackingDetector.h", 19,
                  typeid(::TriNonTrackingDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriNonTrackingDetector::Dictionary, isa_proxy, 4,
                  sizeof(::TriNonTrackingDetector) );
      instance.SetDelete(&delete_TriNonTrackingDetector);
      instance.SetDeleteArray(&deleteArray_TriNonTrackingDetector);
      instance.SetDestructor(&destruct_TriNonTrackingDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriNonTrackingDetector*)
   {
      return GenerateInitInstanceLocal((::TriNonTrackingDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriNonTrackingDetector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriScintillator(void *p = 0);
   static void *newArray_TriScintillator(Long_t size, void *p);
   static void delete_TriScintillator(void *p);
   static void deleteArray_TriScintillator(void *p);
   static void destruct_TriScintillator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriScintillator*)
   {
      ::TriScintillator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriScintillator >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriScintillator", ::TriScintillator::Class_Version(), "TriScintillator.h", 15,
                  typeid(::TriScintillator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriScintillator::Dictionary, isa_proxy, 4,
                  sizeof(::TriScintillator) );
      instance.SetNew(&new_TriScintillator);
      instance.SetNewArray(&newArray_TriScintillator);
      instance.SetDelete(&delete_TriScintillator);
      instance.SetDeleteArray(&deleteArray_TriScintillator);
      instance.SetDestructor(&destruct_TriScintillator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriScintillator*)
   {
      return GenerateInitInstanceLocal((::TriScintillator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriScintillator*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriVDCAnalyticTTDConv(void *p = 0);
   static void *newArray_TriVDCAnalyticTTDConv(Long_t size, void *p);
   static void delete_TriVDCAnalyticTTDConv(void *p);
   static void deleteArray_TriVDCAnalyticTTDConv(void *p);
   static void destruct_TriVDCAnalyticTTDConv(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCAnalyticTTDConv*)
   {
      ::TriVDCAnalyticTTDConv *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCAnalyticTTDConv >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCAnalyticTTDConv", ::TriVDCAnalyticTTDConv::Class_Version(), "TriVDCAnalyticTTDConv.h", 13,
                  typeid(::TriVDCAnalyticTTDConv), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCAnalyticTTDConv::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCAnalyticTTDConv) );
      instance.SetNew(&new_TriVDCAnalyticTTDConv);
      instance.SetNewArray(&newArray_TriVDCAnalyticTTDConv);
      instance.SetDelete(&delete_TriVDCAnalyticTTDConv);
      instance.SetDeleteArray(&deleteArray_TriVDCAnalyticTTDConv);
      instance.SetDestructor(&destruct_TriVDCAnalyticTTDConv);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCAnalyticTTDConv*)
   {
      return GenerateInitInstanceLocal((::TriVDCAnalyticTTDConv*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCAnalyticTTDConv*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriTriggerTime(void *p = 0);
   static void *newArray_TriTriggerTime(Long_t size, void *p);
   static void delete_TriTriggerTime(void *p);
   static void deleteArray_TriTriggerTime(void *p);
   static void destruct_TriTriggerTime(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriTriggerTime*)
   {
      ::TriTriggerTime *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriTriggerTime >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriTriggerTime", ::TriTriggerTime::Class_Version(), "TriTriggerTime.h", 19,
                  typeid(::TriTriggerTime), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriTriggerTime::Dictionary, isa_proxy, 4,
                  sizeof(::TriTriggerTime) );
      instance.SetNew(&new_TriTriggerTime);
      instance.SetNewArray(&newArray_TriTriggerTime);
      instance.SetDelete(&delete_TriTriggerTime);
      instance.SetDeleteArray(&deleteArray_TriTriggerTime);
      instance.SetDestructor(&destruct_TriTriggerTime);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriTriggerTime*)
   {
      return GenerateInitInstanceLocal((::TriTriggerTime*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriTriggerTime*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriTrackProj(void *p = 0);
   static void *newArray_TriTrackProj(Long_t size, void *p);
   static void delete_TriTrackProj(void *p);
   static void deleteArray_TriTrackProj(void *p);
   static void destruct_TriTrackProj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriTrackProj*)
   {
      ::TriTrackProj *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriTrackProj >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriTrackProj", ::TriTrackProj::Class_Version(), "TriTrackProj.h", 16,
                  typeid(::TriTrackProj), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriTrackProj::Dictionary, isa_proxy, 4,
                  sizeof(::TriTrackProj) );
      instance.SetNew(&new_TriTrackProj);
      instance.SetNewArray(&newArray_TriTrackProj);
      instance.SetDelete(&delete_TriTrackProj);
      instance.SetDeleteArray(&deleteArray_TriTrackProj);
      instance.SetDestructor(&destruct_TriTrackProj);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriTrackProj*)
   {
      return GenerateInitInstanceLocal((::TriTrackProj*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriTrackProj*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TriXscin(void *p = 0);
   static void *newArray_TriXscin(Long_t size, void *p);
   static void delete_TriXscin(void *p);
   static void deleteArray_TriXscin(void *p);
   static void destruct_TriXscin(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriXscin*)
   {
      ::TriXscin *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriXscin >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriXscin", ::TriXscin::Class_Version(), "TriXscin.h", 15,
                  typeid(::TriXscin), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriXscin::Dictionary, isa_proxy, 4,
                  sizeof(::TriXscin) );
      instance.SetNew(&new_TriXscin);
      instance.SetNewArray(&newArray_TriXscin);
      instance.SetDelete(&delete_TriXscin);
      instance.SetDeleteArray(&deleteArray_TriXscin);
      instance.SetDestructor(&destruct_TriXscin);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriXscin*)
   {
      return GenerateInitInstanceLocal((::TriXscin*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriXscin*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TriTrackInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriTrackInfo::Class_Name()
{
   return "TriTrackInfo";
}

//______________________________________________________________________________
const char *TriTrackInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriTrackInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriTrackInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriTrackInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriTrackingModule::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriTrackingModule::Class_Name()
{
   return "TriTrackingModule";
}

//______________________________________________________________________________
const char *TriTrackingModule::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingModule*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriTrackingModule::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingModule*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriTrackingModule::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingModule*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriTrackingModule::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingModule*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriSpectrometer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriSpectrometer::Class_Name()
{
   return "TriSpectrometer";
}

//______________________________________________________________________________
const char *TriSpectrometer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriSpectrometer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriSpectrometer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriSpectrometer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriHRS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriHRS::Class_Name()
{
   return "TriHRS";
}

//______________________________________________________________________________
const char *TriHRS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriHRS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriHRS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriHRS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriHRS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriHRS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriHRS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriHRS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriSpectrometerDetector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriSpectrometerDetector::Class_Name()
{
   return "TriSpectrometerDetector";
}

//______________________________________________________________________________
const char *TriSpectrometerDetector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometerDetector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriSpectrometerDetector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometerDetector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriSpectrometerDetector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometerDetector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriSpectrometerDetector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriSpectrometerDetector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriTrackingDetector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriTrackingDetector::Class_Name()
{
   return "TriTrackingDetector";
}

//______________________________________________________________________________
const char *TriTrackingDetector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingDetector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriTrackingDetector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingDetector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriTrackingDetector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingDetector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriTrackingDetector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackingDetector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDC::Class_Name()
{
   return "TriVDC";
}

//______________________________________________________________________________
const char *TriVDC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCPlane::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCPlane::Class_Name()
{
   return "TriVDCPlane";
}

//______________________________________________________________________________
const char *TriVDCPlane::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCPlane*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCPlane::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCPlane*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCPlane::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCPlane*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCPlane::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCPlane*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCUVPlane::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCUVPlane::Class_Name()
{
   return "TriVDCUVPlane";
}

//______________________________________________________________________________
const char *TriVDCUVPlane::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVPlane*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCUVPlane::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVPlane*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCUVPlane::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVPlane*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCUVPlane::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVPlane*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCCluster::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCCluster::Class_Name()
{
   return "TriVDCCluster";
}

//______________________________________________________________________________
const char *TriVDCCluster::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCCluster*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCCluster::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCCluster*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCCluster::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCCluster*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCCluster::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCCluster*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCWire::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCWire::Class_Name()
{
   return "TriVDCWire";
}

//______________________________________________________________________________
const char *TriVDCWire::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCWire*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCWire::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCWire*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCWire::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCWire*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCWire::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCWire*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCHit::Class_Name()
{
   return "TriVDCHit";
}

//______________________________________________________________________________
const char *TriVDCHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCTimeToDistConv::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCTimeToDistConv::Class_Name()
{
   return "TriVDCTimeToDistConv";
}

//______________________________________________________________________________
const char *TriVDCTimeToDistConv::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTimeToDistConv*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCTimeToDistConv::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTimeToDistConv*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCTimeToDistConv::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTimeToDistConv*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCTimeToDistConv::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTimeToDistConv*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCUVTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCUVTrack::Class_Name()
{
   return "TriVDCUVTrack";
}

//______________________________________________________________________________
const char *TriVDCUVTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCUVTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCUVTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCUVTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCUVTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCTrackPair::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCTrackPair::Class_Name()
{
   return "TriVDCTrackPair";
}

//______________________________________________________________________________
const char *TriVDCTrackPair::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackPair*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCTrackPair::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackPair*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCTrackPair::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackPair*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCTrackPair::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackPair*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriTrackID::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriTrackID::Class_Name()
{
   return "TriTrackID";
}

//______________________________________________________________________________
const char *TriTrackID::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackID*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriTrackID::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackID*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriTrackID::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackID*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriTrackID::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackID*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCTrackID::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCTrackID::Class_Name()
{
   return "TriVDCTrackID";
}

//______________________________________________________________________________
const char *TriVDCTrackID::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackID*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCTrackID::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackID*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCTrackID::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackID*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCTrackID::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCTrackID*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriTrack::Class_Name()
{
   return "TriTrack";
}

//______________________________________________________________________________
const char *TriTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriNonTrackingDetector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriNonTrackingDetector::Class_Name()
{
   return "TriNonTrackingDetector";
}

//______________________________________________________________________________
const char *TriNonTrackingDetector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriNonTrackingDetector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriNonTrackingDetector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriNonTrackingDetector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriNonTrackingDetector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriNonTrackingDetector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriNonTrackingDetector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriNonTrackingDetector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriScintillator::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriScintillator::Class_Name()
{
   return "TriScintillator";
}

//______________________________________________________________________________
const char *TriScintillator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriScintillator*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriScintillator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriScintillator*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriScintillator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriScintillator*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriScintillator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriScintillator*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriVDCAnalyticTTDConv::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCAnalyticTTDConv::Class_Name()
{
   return "TriVDCAnalyticTTDConv";
}

//______________________________________________________________________________
const char *TriVDCAnalyticTTDConv::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCAnalyticTTDConv*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCAnalyticTTDConv::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCAnalyticTTDConv*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCAnalyticTTDConv::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCAnalyticTTDConv*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCAnalyticTTDConv::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCAnalyticTTDConv*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriTriggerTime::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriTriggerTime::Class_Name()
{
   return "TriTriggerTime";
}

//______________________________________________________________________________
const char *TriTriggerTime::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTriggerTime*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriTriggerTime::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTriggerTime*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriTriggerTime::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTriggerTime*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriTriggerTime::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTriggerTime*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriTrackProj::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriTrackProj::Class_Name()
{
   return "TriTrackProj";
}

//______________________________________________________________________________
const char *TriTrackProj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackProj*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriTrackProj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriTrackProj*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriTrackProj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackProj*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriTrackProj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriTrackProj*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TriXscin::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriXscin::Class_Name()
{
   return "TriXscin";
}

//______________________________________________________________________________
const char *TriXscin::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriXscin*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriXscin::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriXscin*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriXscin::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriXscin*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriXscin::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriXscin*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TriTrackInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriTrackInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriTrackInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriTrackInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriTrackInfo(void *p) {
      return  p ? new(p) ::TriTrackInfo : new ::TriTrackInfo;
   }
   static void *newArray_TriTrackInfo(Long_t nElements, void *p) {
      return p ? new(p) ::TriTrackInfo[nElements] : new ::TriTrackInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriTrackInfo(void *p) {
      delete ((::TriTrackInfo*)p);
   }
   static void deleteArray_TriTrackInfo(void *p) {
      delete [] ((::TriTrackInfo*)p);
   }
   static void destruct_TriTrackInfo(void *p) {
      typedef ::TriTrackInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriTrackInfo

//______________________________________________________________________________
void TriTrackingModule::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriTrackingModule.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriTrackingModule::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriTrackingModule::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriTrackingModule(void *p) {
      return  p ? new(p) ::TriTrackingModule : new ::TriTrackingModule;
   }
   static void *newArray_TriTrackingModule(Long_t nElements, void *p) {
      return p ? new(p) ::TriTrackingModule[nElements] : new ::TriTrackingModule[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriTrackingModule(void *p) {
      delete ((::TriTrackingModule*)p);
   }
   static void deleteArray_TriTrackingModule(void *p) {
      delete [] ((::TriTrackingModule*)p);
   }
   static void destruct_TriTrackingModule(void *p) {
      typedef ::TriTrackingModule current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriTrackingModule

//______________________________________________________________________________
void TriSpectrometer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriSpectrometer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriSpectrometer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriSpectrometer::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriSpectrometer(void *p) {
      delete ((::TriSpectrometer*)p);
   }
   static void deleteArray_TriSpectrometer(void *p) {
      delete [] ((::TriSpectrometer*)p);
   }
   static void destruct_TriSpectrometer(void *p) {
      typedef ::TriSpectrometer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriSpectrometer

//______________________________________________________________________________
void TriHRS::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriHRS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriHRS::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriHRS::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriHRS(void *p) {
      delete ((::TriHRS*)p);
   }
   static void deleteArray_TriHRS(void *p) {
      delete [] ((::TriHRS*)p);
   }
   static void destruct_TriHRS(void *p) {
      typedef ::TriHRS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriHRS

//______________________________________________________________________________
void TriSpectrometerDetector::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriSpectrometerDetector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriSpectrometerDetector::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriSpectrometerDetector::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriSpectrometerDetector(void *p) {
      delete ((::TriSpectrometerDetector*)p);
   }
   static void deleteArray_TriSpectrometerDetector(void *p) {
      delete [] ((::TriSpectrometerDetector*)p);
   }
   static void destruct_TriSpectrometerDetector(void *p) {
      typedef ::TriSpectrometerDetector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriSpectrometerDetector

//______________________________________________________________________________
void TriTrackingDetector::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriTrackingDetector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriTrackingDetector::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriTrackingDetector::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriTrackingDetector(void *p) {
      delete ((::TriTrackingDetector*)p);
   }
   static void deleteArray_TriTrackingDetector(void *p) {
      delete [] ((::TriTrackingDetector*)p);
   }
   static void destruct_TriTrackingDetector(void *p) {
      typedef ::TriTrackingDetector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriTrackingDetector

//______________________________________________________________________________
void TriVDC::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDC::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDC::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriVDC(void *p) {
      delete ((::TriVDC*)p);
   }
   static void deleteArray_TriVDC(void *p) {
      delete [] ((::TriVDC*)p);
   }
   static void destruct_TriVDC(void *p) {
      typedef ::TriVDC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDC

//______________________________________________________________________________
void TriVDCPlane::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCPlane.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCPlane::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCPlane::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCPlane(void *p) {
      return  p ? new(p) ::TriVDCPlane : new ::TriVDCPlane;
   }
   static void *newArray_TriVDCPlane(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCPlane[nElements] : new ::TriVDCPlane[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCPlane(void *p) {
      delete ((::TriVDCPlane*)p);
   }
   static void deleteArray_TriVDCPlane(void *p) {
      delete [] ((::TriVDCPlane*)p);
   }
   static void destruct_TriVDCPlane(void *p) {
      typedef ::TriVDCPlane current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCPlane

//______________________________________________________________________________
void TriVDCUVPlane::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCUVPlane.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCUVPlane::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCUVPlane::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCUVPlane(void *p) {
      return  p ? new(p) ::TriVDCUVPlane : new ::TriVDCUVPlane;
   }
   static void *newArray_TriVDCUVPlane(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCUVPlane[nElements] : new ::TriVDCUVPlane[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCUVPlane(void *p) {
      delete ((::TriVDCUVPlane*)p);
   }
   static void deleteArray_TriVDCUVPlane(void *p) {
      delete [] ((::TriVDCUVPlane*)p);
   }
   static void destruct_TriVDCUVPlane(void *p) {
      typedef ::TriVDCUVPlane current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCUVPlane

//______________________________________________________________________________
void TriVDCCluster::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCCluster.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCCluster::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCCluster::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCCluster(void *p) {
      return  p ? new(p) ::TriVDCCluster : new ::TriVDCCluster;
   }
   static void *newArray_TriVDCCluster(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCCluster[nElements] : new ::TriVDCCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCCluster(void *p) {
      delete ((::TriVDCCluster*)p);
   }
   static void deleteArray_TriVDCCluster(void *p) {
      delete [] ((::TriVDCCluster*)p);
   }
   static void destruct_TriVDCCluster(void *p) {
      typedef ::TriVDCCluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCCluster

//______________________________________________________________________________
void TriVDCWire::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCWire.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCWire::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCWire::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCWire(void *p) {
      return  p ? new(p) ::TriVDCWire : new ::TriVDCWire;
   }
   static void *newArray_TriVDCWire(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCWire[nElements] : new ::TriVDCWire[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCWire(void *p) {
      delete ((::TriVDCWire*)p);
   }
   static void deleteArray_TriVDCWire(void *p) {
      delete [] ((::TriVDCWire*)p);
   }
   static void destruct_TriVDCWire(void *p) {
      typedef ::TriVDCWire current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCWire

//______________________________________________________________________________
void TriVDCHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCHit(void *p) {
      return  p ? new(p) ::TriVDCHit : new ::TriVDCHit;
   }
   static void *newArray_TriVDCHit(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCHit[nElements] : new ::TriVDCHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCHit(void *p) {
      delete ((::TriVDCHit*)p);
   }
   static void deleteArray_TriVDCHit(void *p) {
      delete [] ((::TriVDCHit*)p);
   }
   static void destruct_TriVDCHit(void *p) {
      typedef ::TriVDCHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCHit

//______________________________________________________________________________
void TriVDCTimeToDistConv::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCTimeToDistConv.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCTimeToDistConv::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCTimeToDistConv::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriVDCTimeToDistConv(void *p) {
      delete ((::TriVDCTimeToDistConv*)p);
   }
   static void deleteArray_TriVDCTimeToDistConv(void *p) {
      delete [] ((::TriVDCTimeToDistConv*)p);
   }
   static void destruct_TriVDCTimeToDistConv(void *p) {
      typedef ::TriVDCTimeToDistConv current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCTimeToDistConv

//______________________________________________________________________________
void TriVDCUVTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCUVTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCUVTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCUVTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCUVTrack(void *p) {
      return  p ? new(p) ::TriVDCUVTrack : new ::TriVDCUVTrack;
   }
   static void *newArray_TriVDCUVTrack(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCUVTrack[nElements] : new ::TriVDCUVTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCUVTrack(void *p) {
      delete ((::TriVDCUVTrack*)p);
   }
   static void deleteArray_TriVDCUVTrack(void *p) {
      delete [] ((::TriVDCUVTrack*)p);
   }
   static void destruct_TriVDCUVTrack(void *p) {
      typedef ::TriVDCUVTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCUVTrack

//______________________________________________________________________________
void TriVDCTrackPair::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCTrackPair.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCTrackPair::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCTrackPair::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCTrackPair(void *p) {
      return  p ? new(p) ::TriVDCTrackPair : new ::TriVDCTrackPair;
   }
   static void *newArray_TriVDCTrackPair(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCTrackPair[nElements] : new ::TriVDCTrackPair[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCTrackPair(void *p) {
      delete ((::TriVDCTrackPair*)p);
   }
   static void deleteArray_TriVDCTrackPair(void *p) {
      delete [] ((::TriVDCTrackPair*)p);
   }
   static void destruct_TriVDCTrackPair(void *p) {
      typedef ::TriVDCTrackPair current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCTrackPair

//______________________________________________________________________________
void TriTrackID::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriTrackID.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriTrackID::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriTrackID::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriTrackID(void *p) {
      delete ((::TriTrackID*)p);
   }
   static void deleteArray_TriTrackID(void *p) {
      delete [] ((::TriTrackID*)p);
   }
   static void destruct_TriTrackID(void *p) {
      typedef ::TriTrackID current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriTrackID

//______________________________________________________________________________
void TriVDCTrackID::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCTrackID.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCTrackID::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCTrackID::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCTrackID(void *p) {
      return  p ? new(p) ::TriVDCTrackID : new ::TriVDCTrackID;
   }
   static void *newArray_TriVDCTrackID(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCTrackID[nElements] : new ::TriVDCTrackID[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCTrackID(void *p) {
      delete ((::TriVDCTrackID*)p);
   }
   static void deleteArray_TriVDCTrackID(void *p) {
      delete [] ((::TriVDCTrackID*)p);
   }
   static void destruct_TriVDCTrackID(void *p) {
      typedef ::TriVDCTrackID current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCTrackID

//______________________________________________________________________________
void TriTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriTrack(void *p) {
      return  p ? new(p) ::TriTrack : new ::TriTrack;
   }
   static void *newArray_TriTrack(Long_t nElements, void *p) {
      return p ? new(p) ::TriTrack[nElements] : new ::TriTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriTrack(void *p) {
      delete ((::TriTrack*)p);
   }
   static void deleteArray_TriTrack(void *p) {
      delete [] ((::TriTrack*)p);
   }
   static void destruct_TriTrack(void *p) {
      typedef ::TriTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriTrack

//______________________________________________________________________________
void TriNonTrackingDetector::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriNonTrackingDetector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriNonTrackingDetector::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriNonTrackingDetector::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriNonTrackingDetector(void *p) {
      delete ((::TriNonTrackingDetector*)p);
   }
   static void deleteArray_TriNonTrackingDetector(void *p) {
      delete [] ((::TriNonTrackingDetector*)p);
   }
   static void destruct_TriNonTrackingDetector(void *p) {
      typedef ::TriNonTrackingDetector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriNonTrackingDetector

//______________________________________________________________________________
void TriScintillator::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriScintillator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriScintillator::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriScintillator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriScintillator(void *p) {
      return  p ? new(p) ::TriScintillator : new ::TriScintillator;
   }
   static void *newArray_TriScintillator(Long_t nElements, void *p) {
      return p ? new(p) ::TriScintillator[nElements] : new ::TriScintillator[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriScintillator(void *p) {
      delete ((::TriScintillator*)p);
   }
   static void deleteArray_TriScintillator(void *p) {
      delete [] ((::TriScintillator*)p);
   }
   static void destruct_TriScintillator(void *p) {
      typedef ::TriScintillator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriScintillator

//______________________________________________________________________________
void TriVDCAnalyticTTDConv::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCAnalyticTTDConv.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCAnalyticTTDConv::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCAnalyticTTDConv::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriVDCAnalyticTTDConv(void *p) {
      return  p ? new(p) ::TriVDCAnalyticTTDConv : new ::TriVDCAnalyticTTDConv;
   }
   static void *newArray_TriVDCAnalyticTTDConv(Long_t nElements, void *p) {
      return p ? new(p) ::TriVDCAnalyticTTDConv[nElements] : new ::TriVDCAnalyticTTDConv[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriVDCAnalyticTTDConv(void *p) {
      delete ((::TriVDCAnalyticTTDConv*)p);
   }
   static void deleteArray_TriVDCAnalyticTTDConv(void *p) {
      delete [] ((::TriVDCAnalyticTTDConv*)p);
   }
   static void destruct_TriVDCAnalyticTTDConv(void *p) {
      typedef ::TriVDCAnalyticTTDConv current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCAnalyticTTDConv

//______________________________________________________________________________
void TriTriggerTime::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriTriggerTime.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriTriggerTime::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriTriggerTime::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriTriggerTime(void *p) {
      return  p ? new(p) ::TriTriggerTime : new ::TriTriggerTime;
   }
   static void *newArray_TriTriggerTime(Long_t nElements, void *p) {
      return p ? new(p) ::TriTriggerTime[nElements] : new ::TriTriggerTime[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriTriggerTime(void *p) {
      delete ((::TriTriggerTime*)p);
   }
   static void deleteArray_TriTriggerTime(void *p) {
      delete [] ((::TriTriggerTime*)p);
   }
   static void destruct_TriTriggerTime(void *p) {
      typedef ::TriTriggerTime current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriTriggerTime

//______________________________________________________________________________
void TriTrackProj::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriTrackProj.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriTrackProj::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriTrackProj::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriTrackProj(void *p) {
      return  p ? new(p) ::TriTrackProj : new ::TriTrackProj;
   }
   static void *newArray_TriTrackProj(Long_t nElements, void *p) {
      return p ? new(p) ::TriTrackProj[nElements] : new ::TriTrackProj[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriTrackProj(void *p) {
      delete ((::TriTrackProj*)p);
   }
   static void deleteArray_TriTrackProj(void *p) {
      delete [] ((::TriTrackProj*)p);
   }
   static void destruct_TriTrackProj(void *p) {
      typedef ::TriTrackProj current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriTrackProj

//______________________________________________________________________________
void TriXscin::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriXscin.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriXscin::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriXscin::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TriXscin(void *p) {
      return  p ? new(p) ::TriXscin : new ::TriXscin;
   }
   static void *newArray_TriXscin(Long_t nElements, void *p) {
      return p ? new(p) ::TriXscin[nElements] : new ::TriXscin[nElements];
   }
   // Wrapper around operator delete
   static void delete_TriXscin(void *p) {
      delete ((::TriXscin*)p);
   }
   static void deleteArray_TriXscin(void *p) {
      delete [] ((::TriXscin*)p);
   }
   static void destruct_TriXscin(void *p) {
      typedef ::TriXscin current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriXscin

namespace {
  void TriggerDictionaryInitialization_TriOldTrackDict_Impl() {
    static const char* headers[] = {
"TriHRS.h",
"TriSpectrometer.h",
"TriSpectrometerDetector.h",
"TriTrackingDetector.h",
"TriVDC.h",
"TriVDCUVPlane.h",
"TriVDCPlane.h",
"TriVDCCluster.h",
"TriVDCHit.h",
"TriVDCTimeToDistConv.h",
"TriVDCUVTrack.h",
"TriVDCTrackPair.h",
"TriVDCTrackID.h",
"TriVDCWire.h",
"TriTrack.h",
"TriTrackInfo.h",
"TriTrackingModule.h",
"TriNonTrackingDetector.h",
"TriScintillator.h",
"TriTrackID.h",
"TriVDCAnalyticTTDConv.h",
"TriTriggerTime.h",
"TriTrackProj.h",
"TriXscin.h",
0
    };
    static const char* includePaths[] = {
"/u/apps/root/6.10.02/root/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/src",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/hana_decode",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriOldTrack",
"/u/apps/root/6.10.02/root/include",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriOldTrack/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TriOldTrackDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Track information)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriTrackInfo.h")))  __attribute__((annotate("$clingAutoload$TriHRS.h")))  TriTrackInfo;
class __attribute__((annotate(R"ATTRDUMP(ABC for a tracking module)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriTrackingModule.h")))  __attribute__((annotate("$clingAutoload$TriHRS.h")))  TriTrackingModule;
class __attribute__((annotate(R"ATTRDUMP(A generic spectrometer)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriSpectrometer.h")))  __attribute__((annotate("$clingAutoload$TriHRS.h")))  TriSpectrometer;
class __attribute__((annotate(R"ATTRDUMP(A Hall A High Resolution Spectrometer)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriHRS.h")))  TriHRS;
class __attribute__((annotate(R"ATTRDUMP(ABC for a spectrometer detector)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriSpectrometerDetector.h")))  TriSpectrometerDetector;
class __attribute__((annotate(R"ATTRDUMP(ABC for a generic tracking detector)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriTrackingDetector.h")))  TriTrackingDetector;
class __attribute__((annotate(R"ATTRDUMP(VDC class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDC.h")))  TriVDC;
class __attribute__((annotate(R"ATTRDUMP(VDCPlane class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCPlane.h")))  __attribute__((annotate("$clingAutoload$TriVDCUVPlane.h")))  TriVDCPlane;
class __attribute__((annotate(R"ATTRDUMP(VDCUVPlane class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCUVPlane.h")))  TriVDCUVPlane;
class __attribute__((annotate(R"ATTRDUMP(A group of VDC hits)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCCluster.h")))  TriVDCCluster;
class __attribute__((annotate(R"ATTRDUMP(VDCWire class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCWire.h")))  __attribute__((annotate("$clingAutoload$TriVDCHit.h")))  TriVDCWire;
class __attribute__((annotate(R"ATTRDUMP(VDCHit class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCHit.h")))  TriVDCHit;
class __attribute__((annotate(R"ATTRDUMP(VDCTimeToDistConv class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCTimeToDistConv.h")))  TriVDCTimeToDistConv;
class __attribute__((annotate(R"ATTRDUMP(VDCUVTrack class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCUVTrack.h")))  TriVDCUVTrack;
class __attribute__((annotate(R"ATTRDUMP(A pair of VDC UV tracks)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCTrackPair.h")))  TriVDCTrackPair;
class __attribute__((annotate(R"ATTRDUMP(Track ID abstract base class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriTrackID.h")))  __attribute__((annotate("$clingAutoload$TriVDCTrackID.h")))  TriTrackID;
class __attribute__((annotate(R"ATTRDUMP(Track ID class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCTrackID.h")))  TriVDCTrackID;
class __attribute__((annotate(R"ATTRDUMP(A generic particle track)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriTrack.h")))  TriTrack;
class __attribute__((annotate(R"ATTRDUMP(ABC for a non-tracking spectrometer detector)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriNonTrackingDetector.h")))  TriNonTrackingDetector;
class __attribute__((annotate(R"ATTRDUMP(Generic scintillator class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriScintillator.h")))  TriScintillator;
class __attribute__((annotate(R"ATTRDUMP(VDC Analytic TTD Conv class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCAnalyticTTDConv.h")))  TriVDCAnalyticTTDConv;
class __attribute__((annotate("$clingAutoload$TriTriggerTime.h")))  TriTriggerTime;
class __attribute__((annotate(R"ATTRDUMP(track projected to detector plane)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriTrackProj.h")))  TriTrackProj;
class __attribute__((annotate(R"ATTRDUMP(Generic scintillator class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriXscin.h")))  TriXscin;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TriOldTrackDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef NDEBUG
  #define NDEBUG 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TriHRS.h"
#include "TriSpectrometer.h"
#include "TriSpectrometerDetector.h"
#include "TriTrackingDetector.h"
#include "TriVDC.h"
#include "TriVDCUVPlane.h"
#include "TriVDCPlane.h"
#include "TriVDCCluster.h"
#include "TriVDCHit.h"
#include "TriVDCTimeToDistConv.h"
#include "TriVDCUVTrack.h"
#include "TriVDCTrackPair.h"
#include "TriVDCTrackID.h"
#include "TriVDCWire.h"
#include "TriTrack.h"
#include "TriTrackInfo.h"
#include "TriTrackingModule.h"
#include "TriNonTrackingDetector.h"
#include "TriScintillator.h"
#include "TriTrackID.h"
#include "TriVDCAnalyticTTDConv.h"
#include "TriTriggerTime.h"
#include "TriTrackProj.h"
#include "TriXscin.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TriHRS", payloadCode, "@",
"TriNonTrackingDetector", payloadCode, "@",
"TriScintillator", payloadCode, "@",
"TriSpectrometer", payloadCode, "@",
"TriSpectrometerDetector", payloadCode, "@",
"TriTrack", payloadCode, "@",
"TriTrackID", payloadCode, "@",
"TriTrackInfo", payloadCode, "@",
"TriTrackProj", payloadCode, "@",
"TriTrackingDetector", payloadCode, "@",
"TriTrackingModule", payloadCode, "@",
"TriTriggerTime", payloadCode, "@",
"TriVDC", payloadCode, "@",
"TriVDCAnalyticTTDConv", payloadCode, "@",
"TriVDCCluster", payloadCode, "@",
"TriVDCHit", payloadCode, "@",
"TriVDCPlane", payloadCode, "@",
"TriVDCTimeToDistConv", payloadCode, "@",
"TriVDCTrackID", payloadCode, "@",
"TriVDCTrackPair", payloadCode, "@",
"TriVDCUVPlane", payloadCode, "@",
"TriVDCUVTrack", payloadCode, "@",
"TriVDCWire", payloadCode, "@",
"TriXscin", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TriOldTrackDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TriOldTrackDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TriOldTrackDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TriOldTrackDict() {
  TriggerDictionaryInitialization_TriOldTrackDict_Impl();
}
