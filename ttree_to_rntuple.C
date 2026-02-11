// Heavily inspired from this tutorial: https://root.cern/doc/v636/ntpl008__import_8C.html

#include <ROOT/RNTupleDS.hxx>
#include <ROOT/RNTupleImporter.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RPageStorageFile.hxx>
 
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
 
using RNTupleImporter = ROOT::Experimental::RNTupleImporter;

constexpr char const *kTreeName = "FlatTree_VARS";
constexpr char const *kRNTupleName = "Events";

int main(int argc, char const *argv[])
{
  if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <tree-file-name>  <ntuple-file-name>" << std::endl;
      return 1;
  }
  char const *kTreeFileName = argv[1];
  char const *kNTupleFileName = argv[2];

  // RNTupleImporter appends keys to the output file; make sure a second run of the tutorial does not fail
  // with `Key 'Events' already exists in file ntpl008_import.root` by removing the output file.
  gSystem->Unlink(kNTupleFileName);
 
  ROOT::EnableImplicitMT();
  auto importer = RNTupleImporter::Create(kTreeFileName, kTreeName, kNTupleFileName);
  importer->SetNTupleName(kRNTupleName);
  importer->Import();
  return 0;
}