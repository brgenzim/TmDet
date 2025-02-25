<?php

namespace Tmdet\Compare\Xml;

use SimpleXMLElement;
use Tmdet\Compare\Helper;
use Tmdet\Compare\Xml\Constant4 as XML;

class Xml4 {
    private $xml;
    private $path;

    public function __construct(string $path) {
        if (file_exists($path)) {
            $this->xml = simplexml_load_file($path);
        }
        else {
            $this->xml = false;
        }
        $this->path = $path;
    }

    public function getTmpStatus() {
        return $this->xml ? $this->xml[XML::ATTR_TRANSMEMBRANE]->__toString() : "missing";
    }

    public function getProteinType() {
        return $this->xml ? $this->xml->{XML::NODE_RAWDATA}->{ Xml::NODE_TMTYPE }->__toString() : "missing";
    }

    public function getChains() : array {
        $ret = [];
        if ($this->xml !== false 
            && $this->xml !== null 
            && $this->xml->{XML::NODE_CHAINS} !== null
            && $this->xml->{XML::NODE_CHAINS}->{XML::NODE_CHAIN} !== null) {
            foreach($this->xml->{XML::NODE_CHAINS}->{XML::NODE_CHAIN} as $chain) {
                $id = $chain[XML::ATTR_AUTH_ID]->__toString();
                $id = ($id=="_"?"A":$id);
                $type = $chain[XML::ATTR_TYPE]->__toString();
                //if ($type == 'alpha' || $type == 'beta') {
                    $ret[$id] = [
                        'id' => $id,
                        'type' => $type,
                        'numtm' => $chain[XML::ATTR_NUM_TM]->__toString(),
                        'seq' => Helper::clearSeq($chain->{XML::NODE_SEQENCE}->__toString()),
                        'regions' => $this->getRegions($chain),
                    ];
                //}
            }
        }
        else {
            echo "Problem with the xml file: ".$this->path."\n";
        }
        return $ret;
    }

    private function getRegions(SimpleXMLElement $chainNode) : array {
        $ret = [];
        if ( $chainNode !== null
            && $chainNode->{XML::NODE_REGIONS} !== null
            && $chainNode->{XML::NODE_REGIONS}->{XML::NODE_REGION} !== null) {
            foreach($chainNode->{XML::NODE_REGIONS}->{XML::NODE_REGION} as $region) {
                $ret[] = [
                    'start' => intval($region[XML::ATTR_START_LABEL_ID])-1,
                    'end' => intval($region[XML::ATTR_END_LABEL_ID])-1,
                    'type' => $region[XML::ATTR_TYPE]->__toString(),
                    'astart' => intval($region[XML::ATTR_START_AUTH_ID])-1,
                    'aend' => intval($region[XML::ATTR_END_AUTH_ID])-1,
                ];
            }
        }
        return $ret;
    }

    public function hasRequestedRegion(string $regionCode) : bool {
        foreach($this->getChains() as $chain) {
            foreach($chain['regions'] as $region) {
                if ($region['type'] == $regionCode) {
                    return true;
                }
            }
        }
        return false;
    }
}
