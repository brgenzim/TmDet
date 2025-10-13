// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <Config.hpp>
#include <DTOs/Xml.hpp>
#include <DTOs/XmlRW/Reader.hpp>
#include <DTOs/XmlRW/Writer.hpp>
#include <DTOs/XmlRW/Writer3.hpp>
#include <System/Arguments.hpp>
#include <System/Logger.hpp>
#include <System/FilePaths.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Xml.hpp>

namespace Tmdet::DTOs {

        void Xml::fromProtein(const Tmdet::VOs::Protein& protein) {
            xmlData.code = protein.code;
            xmlData.tmp = protein.tmp;
            xmlData.date = protein.date;
            xmlData.version = protein.version;
            xmlData.modifications = protein.modifications;
            xmlData.qValue = protein.qValue;
            xmlData.type = protein.type;
            xmlData.spres = protein.spres;
            xmlData.pdbkwres = protein.pdbkwres;
            xmlData.bioMatrix = protein.bioMatrix;
            xmlData.membranes = protein.membranes;
            xmlData.tmatrix = protein.tmatrix;
            for(const auto& chain: protein.chains) {
                if (chain.labelId != "") {
                    xmlData.chains.emplace_back(chain.id,
                        chain.labelId,
                        chain.selected,
                        chain.numtm,
                        chain.seq,
                        chain.regions,
                        chain.type);
                }
            }
        }

        void Xml::toProtein(Tmdet::VOs::Protein& protein) const {
            protein.code = xmlData.code;
            protein.tmp = xmlData.tmp;
            protein.date = xmlData.date;
            protein.version = xmlData.version;
            protein.modifications = xmlData.modifications;
            protein.qValue = xmlData.qValue;
            protein.type = xmlData.type;
            protein.spres = xmlData.spres;
            protein.pdbkwres = xmlData.pdbkwres;
            protein.bioMatrix = xmlData.bioMatrix;
            protein.membranes = xmlData.membranes;
            protein.tmatrix = xmlData.tmatrix;
            for(const auto& xmlChain: xmlData.chains) {
                bool found = false;
                for( auto& proteinChain: protein.chains) {
                    if (xmlChain.id == proteinChain.id) {
                        proteinChain.selected = xmlChain.selected;
                        proteinChain.numtm = xmlChain.numtm;
                        proteinChain.seq = xmlChain.seq; //it should be checked later
                        proteinChain.type = xmlChain.type;
                        proteinChain.regions = xmlChain.regions;
                        found = true;
                        continue;
                    }
                }
                if (!found) {
                    //T O D O: in case of pdbtm 3.0 xml new chain created in biomatrix
                    WARN_LOG("Chain {} not found in pdb structure",xmlChain.id);
                }
            }
        }

        bool Xml::read(const std::string& xmlPath) {
            if ( !Tmdet::System::FilePaths::fileExists(xmlPath) ) {
                WARN_LOG("file not found: {}",xmlPath);
                return false;
            }
            else {
                Tmdet::DTOs::XmlRW::Reader reader;
                reader.readXml(xmlData, xmlPath);
            }
            return true;
        }

        void Xml::read(const std::string& xmlPath, Tmdet::VOs::Protein& protein) {
            if (read(xmlPath)) {
                toProtein(protein);
            }
        }

        void Xml::write(const std::string& xmlPath, const Tmdet::System::Arguments& args) {
            if (outXmlFmt == "v4") {
                Tmdet::DTOs::XmlRW::Writer writer;
                writer.writeXml(xmlData, xmlPath, args);
            }
            else {
                Tmdet::DTOs::XmlRW::Writer3 writer;
                writer.writeXml(xmlData, xmlPath, args);
            }
        }

        void Xml::write(const std::string& xmlPath, const Tmdet::VOs::Protein& protein, const Tmdet::System::Arguments& args) {
            fromProtein(protein);
            write(xmlPath, args);
        }

        void Xml::notTransmembrane(const std::string& xmlInputPath, const std::string& xmlOutputPath, const Tmdet::System::Arguments& args) {
            read(xmlInputPath);
            if (xmlData.version != "") {
                xmlData.notTransmembrane();
                write(xmlOutputPath, args);
            }
            else {
                ERROR_LOG("Cannot overwrite old (version < 4.0) xml file");
            }
        }

        std::string Xml::setPath(const std::string& code, const std::string& x1, const std::string& x2) const {
            if (code != "") {
                return Tmdet::System::FilePaths::xml(code,true);
            }
            return (x1==""?x2:x1);
        }

        void Xml::setV3Fmt() {
            outXmlFmt = "v3";
        }
}
