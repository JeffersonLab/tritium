// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TriVDCeffDict

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
#include "TriVDCeff.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_TriVDCeff(void *p);
   static void deleteArray_TriVDCeff(void *p);
   static void destruct_TriVDCeff(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriVDCeff*)
   {
      ::TriVDCeff *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriVDCeff >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriVDCeff", ::TriVDCeff::Class_Version(), "TriVDCeff.h", 16,
                  typeid(::TriVDCeff), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriVDCeff::Dictionary, isa_proxy, 4,
                  sizeof(::TriVDCeff) );
      instance.SetDelete(&delete_TriVDCeff);
      instance.SetDeleteArray(&deleteArray_TriVDCeff);
      instance.SetDestructor(&destruct_TriVDCeff);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriVDCeff*)
   {
      return GenerateInitInstanceLocal((::TriVDCeff*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriVDCeff*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TriVDCeff::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriVDCeff::Class_Name()
{
   return "TriVDCeff";
}

//______________________________________________________________________________
const char *TriVDCeff::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCeff*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriVDCeff::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriVDCeff*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriVDCeff::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCeff*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriVDCeff::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriVDCeff*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TriVDCeff::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriVDCeff.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriVDCeff::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriVDCeff::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriVDCeff(void *p) {
      delete ((::TriVDCeff*)p);
   }
   static void deleteArray_TriVDCeff(void *p) {
      delete [] ((::TriVDCeff*)p);
   }
   static void destruct_TriVDCeff(void *p) {
      typedef ::TriVDCeff current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriVDCeff

namespace {
  void TriggerDictionaryInitialization_TriVDCeffDict_Impl() {
    static const char* headers[] = {
"TriVDCeff.h",
0
    };
    static const char* includePaths[] = {
"/u/apps/root/6.10.02/root/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/src",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/hana_decode",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriVDCeff",
"/u/apps/root/6.10.02/root/include",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriVDCeff/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TriVDCeffDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(VDC hit efficiency physics module)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriVDCeff.h")))  TriVDCeff;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TriVDCeffDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TriVDCeff.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TriVDCeff", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TriVDCeffDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TriVDCeffDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TriVDCeffDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TriVDCeffDict() {
  TriggerDictionaryInitialization_TriVDCeffDict_Impl();
}
