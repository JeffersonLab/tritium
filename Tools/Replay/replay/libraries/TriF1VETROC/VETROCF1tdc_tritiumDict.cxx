// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME VETROCF1tdc_tritiumDict

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
#include "VETROCtdcModule.h"
#include "TstF1TDCModule.h"
#include "TdcDataEvtHandler.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_DecodercLcLVETROCtdcModule(void *p = 0);
   static void *newArray_DecodercLcLVETROCtdcModule(Long_t size, void *p);
   static void delete_DecodercLcLVETROCtdcModule(void *p);
   static void deleteArray_DecodercLcLVETROCtdcModule(void *p);
   static void destruct_DecodercLcLVETROCtdcModule(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Decoder::VETROCtdcModule*)
   {
      ::Decoder::VETROCtdcModule *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Decoder::VETROCtdcModule >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Decoder::VETROCtdcModule", ::Decoder::VETROCtdcModule::Class_Version(), "VETROCtdcModule.h", 18,
                  typeid(::Decoder::VETROCtdcModule), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Decoder::VETROCtdcModule::Dictionary, isa_proxy, 4,
                  sizeof(::Decoder::VETROCtdcModule) );
      instance.SetNew(&new_DecodercLcLVETROCtdcModule);
      instance.SetNewArray(&newArray_DecodercLcLVETROCtdcModule);
      instance.SetDelete(&delete_DecodercLcLVETROCtdcModule);
      instance.SetDeleteArray(&deleteArray_DecodercLcLVETROCtdcModule);
      instance.SetDestructor(&destruct_DecodercLcLVETROCtdcModule);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Decoder::VETROCtdcModule*)
   {
      return GenerateInitInstanceLocal((::Decoder::VETROCtdcModule*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Decoder::VETROCtdcModule*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_DecodercLcLTstF1TDCModule(void *p = 0);
   static void *newArray_DecodercLcLTstF1TDCModule(Long_t size, void *p);
   static void delete_DecodercLcLTstF1TDCModule(void *p);
   static void deleteArray_DecodercLcLTstF1TDCModule(void *p);
   static void destruct_DecodercLcLTstF1TDCModule(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Decoder::TstF1TDCModule*)
   {
      ::Decoder::TstF1TDCModule *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Decoder::TstF1TDCModule >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Decoder::TstF1TDCModule", ::Decoder::TstF1TDCModule::Class_Version(), "TstF1TDCModule.h", 15,
                  typeid(::Decoder::TstF1TDCModule), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Decoder::TstF1TDCModule::Dictionary, isa_proxy, 4,
                  sizeof(::Decoder::TstF1TDCModule) );
      instance.SetNew(&new_DecodercLcLTstF1TDCModule);
      instance.SetNewArray(&newArray_DecodercLcLTstF1TDCModule);
      instance.SetDelete(&delete_DecodercLcLTstF1TDCModule);
      instance.SetDeleteArray(&deleteArray_DecodercLcLTstF1TDCModule);
      instance.SetDestructor(&destruct_DecodercLcLTstF1TDCModule);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Decoder::TstF1TDCModule*)
   {
      return GenerateInitInstanceLocal((::Decoder::TstF1TDCModule*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Decoder::TstF1TDCModule*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TdcDataEvtHandler(void *p);
   static void deleteArray_TdcDataEvtHandler(void *p);
   static void destruct_TdcDataEvtHandler(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TdcDataEvtHandler*)
   {
      ::TdcDataEvtHandler *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TdcDataEvtHandler >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TdcDataEvtHandler", ::TdcDataEvtHandler::Class_Version(), "TdcDataEvtHandler.h", 19,
                  typeid(::TdcDataEvtHandler), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TdcDataEvtHandler::Dictionary, isa_proxy, 4,
                  sizeof(::TdcDataEvtHandler) );
      instance.SetDelete(&delete_TdcDataEvtHandler);
      instance.SetDeleteArray(&deleteArray_TdcDataEvtHandler);
      instance.SetDestructor(&destruct_TdcDataEvtHandler);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TdcDataEvtHandler*)
   {
      return GenerateInitInstanceLocal((::TdcDataEvtHandler*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TdcDataEvtHandler*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace Decoder {
//______________________________________________________________________________
atomic_TClass_ptr VETROCtdcModule::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VETROCtdcModule::Class_Name()
{
   return "Decoder::VETROCtdcModule";
}

//______________________________________________________________________________
const char *VETROCtdcModule::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Decoder::VETROCtdcModule*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VETROCtdcModule::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Decoder::VETROCtdcModule*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VETROCtdcModule::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Decoder::VETROCtdcModule*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VETROCtdcModule::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Decoder::VETROCtdcModule*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace Decoder
namespace Decoder {
//______________________________________________________________________________
atomic_TClass_ptr TstF1TDCModule::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TstF1TDCModule::Class_Name()
{
   return "Decoder::TstF1TDCModule";
}

//______________________________________________________________________________
const char *TstF1TDCModule::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Decoder::TstF1TDCModule*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TstF1TDCModule::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Decoder::TstF1TDCModule*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TstF1TDCModule::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Decoder::TstF1TDCModule*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TstF1TDCModule::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Decoder::TstF1TDCModule*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace Decoder
//______________________________________________________________________________
atomic_TClass_ptr TdcDataEvtHandler::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TdcDataEvtHandler::Class_Name()
{
   return "TdcDataEvtHandler";
}

//______________________________________________________________________________
const char *TdcDataEvtHandler::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TdcDataEvtHandler*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TdcDataEvtHandler::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TdcDataEvtHandler*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TdcDataEvtHandler::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TdcDataEvtHandler*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TdcDataEvtHandler::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TdcDataEvtHandler*)0x0)->GetClass(); }
   return fgIsA;
}

namespace Decoder {
//______________________________________________________________________________
void VETROCtdcModule::Streamer(TBuffer &R__b)
{
   // Stream an object of class Decoder::VETROCtdcModule.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Decoder::VETROCtdcModule::Class(),this);
   } else {
      R__b.WriteClassBuffer(Decoder::VETROCtdcModule::Class(),this);
   }
}

} // namespace Decoder
namespace ROOT {
   // Wrappers around operator new
   static void *new_DecodercLcLVETROCtdcModule(void *p) {
      return  p ? new(p) ::Decoder::VETROCtdcModule : new ::Decoder::VETROCtdcModule;
   }
   static void *newArray_DecodercLcLVETROCtdcModule(Long_t nElements, void *p) {
      return p ? new(p) ::Decoder::VETROCtdcModule[nElements] : new ::Decoder::VETROCtdcModule[nElements];
   }
   // Wrapper around operator delete
   static void delete_DecodercLcLVETROCtdcModule(void *p) {
      delete ((::Decoder::VETROCtdcModule*)p);
   }
   static void deleteArray_DecodercLcLVETROCtdcModule(void *p) {
      delete [] ((::Decoder::VETROCtdcModule*)p);
   }
   static void destruct_DecodercLcLVETROCtdcModule(void *p) {
      typedef ::Decoder::VETROCtdcModule current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Decoder::VETROCtdcModule

namespace Decoder {
//______________________________________________________________________________
void TstF1TDCModule::Streamer(TBuffer &R__b)
{
   // Stream an object of class Decoder::TstF1TDCModule.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Decoder::TstF1TDCModule::Class(),this);
   } else {
      R__b.WriteClassBuffer(Decoder::TstF1TDCModule::Class(),this);
   }
}

} // namespace Decoder
namespace ROOT {
   // Wrappers around operator new
   static void *new_DecodercLcLTstF1TDCModule(void *p) {
      return  p ? new(p) ::Decoder::TstF1TDCModule : new ::Decoder::TstF1TDCModule;
   }
   static void *newArray_DecodercLcLTstF1TDCModule(Long_t nElements, void *p) {
      return p ? new(p) ::Decoder::TstF1TDCModule[nElements] : new ::Decoder::TstF1TDCModule[nElements];
   }
   // Wrapper around operator delete
   static void delete_DecodercLcLTstF1TDCModule(void *p) {
      delete ((::Decoder::TstF1TDCModule*)p);
   }
   static void deleteArray_DecodercLcLTstF1TDCModule(void *p) {
      delete [] ((::Decoder::TstF1TDCModule*)p);
   }
   static void destruct_DecodercLcLTstF1TDCModule(void *p) {
      typedef ::Decoder::TstF1TDCModule current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Decoder::TstF1TDCModule

//______________________________________________________________________________
void TdcDataEvtHandler::Streamer(TBuffer &R__b)
{
   // Stream an object of class TdcDataEvtHandler.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TdcDataEvtHandler::Class(),this);
   } else {
      R__b.WriteClassBuffer(TdcDataEvtHandler::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TdcDataEvtHandler(void *p) {
      delete ((::TdcDataEvtHandler*)p);
   }
   static void deleteArray_TdcDataEvtHandler(void *p) {
      delete [] ((::TdcDataEvtHandler*)p);
   }
   static void destruct_TdcDataEvtHandler(void *p) {
      typedef ::TdcDataEvtHandler current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TdcDataEvtHandler

namespace {
  void TriggerDictionaryInitialization_VETROCF1tdc_tritiumDict_Impl() {
    static const char* headers[] = {
"VETROCtdcModule.h",
"TstF1TDCModule.h",
"TdcDataEvtHandler.h",
0
    };
    static const char* includePaths[] = {
"/u/apps/root/6.10.02/root/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/include",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/src",
"/work/halla/triton/MARATHON_Replay/analyzer1.6/hana_decode",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriF1VETROC",
"/u/apps/root/6.10.02/root/include",
"/lustre/expphy/work/halla/triton/MARATHON_Replay/replay/libraries/TriF1VETROC/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "VETROCF1tdc_tritiumDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace Decoder{class __attribute__((annotate(R"ATTRDUMP(JLab F1 TDC Module, test version)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$VETROCtdcModule.h")))  VETROCtdcModule;}
namespace Decoder{class __attribute__((annotate(R"ATTRDUMP(JLab F1 TDC Module, test version)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TstF1TDCModule.h")))  TstF1TDCModule;}
class __attribute__((annotate(R"ATTRDUMP(Analyze some TDC data)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TdcDataEvtHandler.h")))  TdcDataEvtHandler;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "VETROCF1tdc_tritiumDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef NDEBUG
  #define NDEBUG 1
#endif
#ifndef WITH_DEBUG
  #define WITH_DEBUG 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "VETROCtdcModule.h"
#include "TstF1TDCModule.h"
#include "TdcDataEvtHandler.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"Decoder::TstF1TDCModule", payloadCode, "@",
"Decoder::VETROCtdcModule", payloadCode, "@",
"TdcDataEvtHandler", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("VETROCF1tdc_tritiumDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_VETROCF1tdc_tritiumDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_VETROCF1tdc_tritiumDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_VETROCF1tdc_tritiumDict() {
  TriggerDictionaryInitialization_VETROCF1tdc_tritiumDict_Impl();
}
