#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Protein.hpp>

/**
 * @brief namespace Tmdet classes
 * @namespace Tmdet
 */
namespace Tmdet {
    
    /**
     * @brief namespace for data transfer objects
     * @namespace DTOs
     */
    namespace DTOs {

    /**
     * @brief base class for reading and writing xml files
     */
    class Xml {
        private:

            /**
             * @brief pugixml document
             */
            pugi::xml_document _doc;

            /**
             * @brief root node of the xml document
             */
            pugi::xml_node _root;

            /**
             * @brief read transformation matrix from xml document object
             * 
             * @param node 
             * @return TmdetVO::TMatrix 
             */
            TmdetVO::TMatrix getTMatrix(const pugi::xml_node& node) const;

            /**
             * @brief writing transformation matrix into xml document object
             * 
             * @param node 
             * @param tmatrix 
             */
            void setTMatrix(const pugi::xml_node& node, TmdetVO::TMatrix& tmatrix) const;

            /**
             * @brief header part of the xml file
             */
            virtual const std::string _pdbtm_xml=nullptr;
            
            /**
             * @brief membrane part of the xml
             */
            virtual const std::string _membrane_xml=nullptr;

            /**
             * @brief description of a chain
             */
            const std::string _chain_xml=R"(
<CHAIN CHAINID="xxx" NUM_TM="xxx" TYPE="xxx">
  <SEQ>
  </SEQ>
</CHAIN>
)";
            const char* XML_ATTR_CHAINID="CHAINID";
            const char* XML_ATTR_ID="ID";
            const char* XML_ATTR_NEW_CHAINID="NEW_CHAINID";
            const char* XML_ATTR_NUM_TM="NUM_TM";
            const char* XML_ATTR_PDB_BEG="pdb_beg";
            const char* XML_ATTR_PDB_BEGI="pdb_begi";
            const char* XML_ATTR_PDB_END="pdb_end";
            const char* XML_ATTR_PDB_ENDI="pdb_endi";
            const char* XML_ATTR_SEQ_BEG="seq_beg";
            const char* XML_ATTR_SEQ_END="seq_end";
            const char* XML_ATTR_TMP="TMP";
            const char* XML_ATTR_TYPE="TYPE";
            const char* XML_ATTR_type="type";
            const char* XML_ATTR_T="T";
            const char* XML_ATTR_X="X";
            const char* XML_ATTR_Y="Y";
            const char* XML_ATTR_Z="Z";
            
            const char* XML_NODE_APPLY_TO_CHAIN="APPLY_TO_CHAIN";
            const char* XML_NODE_BIOMATRIX="BIOMATRIX";
            const char* XML_NODE_CHAIN="CHAIN";
            const char* XML_NODE_CREATE_DATE="CREATE_DATE";
            const char* XML_NODE_DATE="DATE";
            const char* XML_NODE_DELETE="DELETE";
            const char* XML_NODE_DESCR="DESCR";
            const char* XML_NODE_MATRIX="MATRIX";
            const char* XML_NODE_MEMBRANE="MEMBRANE";
            const char* XML_NODE_MODIFICATION="MODIFICATION";
            const char* XML_NODE_NORMAL="NORMAL";
            const char* XML_NODE_ROOT="pdbtm";
            const char* XML_NODE_PDBKWORD="PDBKWORD";
            const char* XML_NODE_PDBKWRES="PDBKWRES";
            const char* XML_NODE_RAWRES="RAWRES";
            const char* XML_NODE_REGION="REGION";
            const char* XML_NODE_ROWX="ROWX";
            const char* XML_NODE_ROWY="ROWY";
            const char* XML_NODE_ROWZ="ROWZ";
            const char* XML_NODE_SEQ="SEQ";
            const char* XML_NODE_SPRES="SPRES";
            const char* XML_NODE_TMATRIX="TMATRIX";
            const char* XML_NODE_TMRES="TMRES";
            const char* XML_NODE_TMTYPE="TMTYPE";
            const char* XML_NODE_TMDET_VERSION="TMDET_VERSION";

        public:
            
            /**
             * @brief read entire xml file into pugixml doucument object
             * 
             * @param path 
             */
            void read(const std::string& path);

            /**
             * @brief writing pugixml object into file
             * 
             * @param path 
             */
            void write(const std::string& path) const;

            /**
             * @brief create an empty pugixml object
             * 
             */
            virtual void create() = 0;

            /**
             * @brief Get the value of tmp attribute
             * 
             * @return true 
             * @return false 
             */
            bool getTmp() const;

            /**
             * @brief Set the value of tmp attribute
             * 
             * @param tmp 
             */
            void setTmp(const bool& tmp) const;

            /**
             * @brief Get the value of code attribute
             * 
             * @return std::string 
             */
            std::string getCode() const;

            /**
             * @brief Set the value of code attribute
             * 
             * @param code 
             */
            void setCode(const std::string& code) const;

            /**
             * @brief Get the value of CREATE_DATE node
             * 
             * @return std::string 
             */
            std::string getCreateDate() const;

            /**
             * @brief Set the value of CREATE_DATE node
             * 
             * @param date 
             */
            void setCreateDate(const std::string& date) const;

            /**
             * @brief Get the value of TMDET_VERSION node
             * 
             * @return std::string 
             */
            std::string getVersion();

            /**
             * @brief Set the value of TMDET_VERSION node
             * 
             * @param version 
             */
            void setVersion(const std::string& version) const;

            /**
             * @brief Get the content of MODIFICATIONS node
             * 
             * @return std::vector<TmdetVO::Modification> 
             */
            std::vector<TmdetVO::Modification> getModifications() const;

            /**
             * @brief Set the content of MODIFICATIONS node
             * 
             * @param mods 
             */
            void setModifications(const std::vector<TmdetVO::Modification>& mods);

            /**
             * @brief Get the value of Qvalue 
             * 
             * @return double 
             */
            double getQvalue() const;

            /**
             * @brief Set the value of Qvalue 
             * 
             * @param q 
             */
            void setQvalue(const double& q) const;

            /**
             * @brief Get the value of transmembrane type
             * 
             * @return std::string 
             */
            std::string getTmtype() const;

            /**
             * @brief Set the  value of transmembrane type
             * 
             * @param type 
             */
            void setTmtype(const std::string& type) const;

            /**
             * @brief Get the value of SPRES node
             * 
             * @return std::string 
             */
            std::string getSpres() const;

            /**
             * @brief Set the value of SPRES node
             * 
             * @param type 
             */
            void setSpres(const std::string& type) const;

            /**
             * @brief Get the value of PDBKWRES node
             * 
             * @return std::string 
             */
            std::string getPdbkwres() const;

            /**
             * @brief Set the value of PDBKWRES node
             * 
             * @param type 
             */
            void setPdbkwres(const std::string& type) const;

            /**
             * @brief Get the content of BIOMATRIX node
             * 
             * @return TmdetVO::BioMatrix 
             */
            TmdetVO::BioMatrix getBioMatrix() const;

            /**
             * @brief Set the content of BIOMATRIX node
             * 
             * @param bioMatrix 
             */
            void setBioMatrix(TmdetVO::BioMatrix& bioMatrix);

            /**
             * @brief Get the content of MEMBRANE node
             * 
             * @return std::vector<TmdetVO::Membrane> 
             */
            std::vector<TmdetVO::Membrane> getMembranes() const;

            /**
             * @brief Set the content of MEMBRANE node
             * 
             * @param membranes 
             */
            void setMembranes(std::vector<TmdetVO::Membrane>& membranes);

            /**
             * @brief Get the content of CHAIN node
             * 
             * @param chains 
             */
            void getChains(std::vector<TmdetVO::Chain>& chains);

            /**
             * @brief Set the content of CHAIN node
             * 
             * @param chains 
             */
            void setChains(const std::vector<TmdetVO::Chain>& chains);

            /**
             * @brief Get the content of REGION node
             * 
             * @param cnode 
             * @return std::vector<TmdetVO::Region> 
             */
            std::vector<TmdetVO::Region> getRegions(const pugi::xml_node& cnode) const;

            /**
             * @brief Set the content of REGION node
             * 
             * @param pnode 
             * @param regions 
             */
            void setRegions(pugi::xml_node& pnode, const std::vector<TmdetVO::Region>& regions) const;

            /**
             * @brief read tmdet xml file into Tmdet Value Object
             * 
             * @param tmdetVO 
             * @param path 
             */
            void readXml(TmdetVO::Protein& tmdetVO, const std::string& path);

            /**
             * @brief write Tmdet Value Objects out in xml
             * 
             * @param tmdetVO 
             * @param path 
             */
            void writeXml(TmdetVO::Protein& tmdetVO, const std::string& path);
    };
}
