// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TriBCMDict

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
#include "TriBCM.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_TriBCM(void *p);
   static void deleteArray_TriBCM(void *p);
   static void destruct_TriBCM(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TriBCM*)
   {
      ::TriBCM *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TriBCM >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TriBCM", ::TriBCM::Class_Version(), "TriBCM.h", 20,
                  typeid(::TriBCM), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TriBCM::Dictionary, isa_proxy, 4,
                  sizeof(::TriBCM) );
      instance.SetDelete(&delete_TriBCM);
      instance.SetDeleteArray(&deleteArray_TriBCM);
      instance.SetDestructor(&destruct_TriBCM);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TriBCM*)
   {
      return GenerateInitInstanceLocal((::TriBCM*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TriBCM*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TriBCM::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriBCM::Class_Name()
{
   return "TriBCM";
}

//______________________________________________________________________________
const char *TriBCM::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriBCM*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriBCM::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TriBCM*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriBCM::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriBCM*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriBCM::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TriBCM*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TriBCM::Streamer(TBuffer &R__b)
{
   // Stream an object of class TriBCM.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TriBCM::Class(),this);
   } else {
      R__b.WriteClassBuffer(TriBCM::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TriBCM(void *p) {
      delete ((::TriBCM*)p);
   }
   static void deleteArray_TriBCM(void *p) {
      delete [] ((::TriBCM*)p);
   }
   static void destruct_TriBCM(void *p) {
      typedef ::TriBCM current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TriBCM

namespace {
  void TriggerDictionaryInitialization_TriBCMDict_Impl() {
    static const char* headers[] = {
"TriBCM.h",
0
    };
    static const char* includePaths[] = {
"/u/apps/root/6.10.02/root/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/src",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/hana_decode",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriBCM",
"/u/apps/root/6.10.02/root/include",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriBCM/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TriBCMDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Scaler Event handler)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TriBCM.h")))  TriBCM;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TriBCMDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TriBCM.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TriBCM", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TriBCMDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TriBCMDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TriBCMDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TriBCMDict() {
  TriggerDictionaryInitialization_TriBCMDict_Impl();
}
