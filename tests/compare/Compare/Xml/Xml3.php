<?php

namespace Tmdet\Compare\Xml;

use Tmdet\Compare\Helper;
use Tmdet\Compare\Xml\Constant3 as XML;

class Xml3 {
    private $xml;

    public function __construct(string $path) {
        if (file_exists($path)) {
            $this->xml = simplexml_load_file($path);
        }
        else {
            $this->xml = false;
        }
    }

    public function getTmpStatus() {
        return $this->xml ? $this->xml[XML::ATTR_TMP]->__toString() : "missing";
    }

    public function getProteinType() {
        return $this->xml ? $this->xml->{XML::NODE_RAWRES}->{ Xml::NODE_TMTYPE }->__toString() : "missing";
    }

    public function getChains() : array {
        $ret = [];
        if ($this->xml !== false) {
            foreach($this->xml->{XML::NODE_CHAIN} as $chain) {
                $id = $chain[XML::ATTR_CHAINID]->__toString();
                $id = ($id=="_"?"A":$id);
                $type = $chain[XML::ATTR_TYPE]->__toString();
                if ($type == 'alpha' || $type == 'beta') {
                    $ret[$id] = [
                        'id' => $id,
                        'type' => $type,
                        'numtm' => $chain[XML::ATTR_NUM_TM]->__toString(),
                        'seq' => Helper::clearSeq($chain->{XML::NODE_SEQ}->__toString()),
                        'regions' => $this->getRegions($chain),
                    ];
                }
            }
        }
        return $ret;
    }

    private function getRegions($chain) : array {
        $ret = [];
        foreach($chain->{XML::NODE_REGION} as $region) {
            $type = $region[XML::ATTR_type]->__toString();
            $type = ($type=='M'?'H':$type);
            $ret [] = [
                'start' => intval($region[XML::ATTR_SEQ_BEG])-1,
                'end' => intval($region[XML::ATTR_SEQ_END])-1,
                'type' => $type,
                'astart' => intval($region[XML::ATTR_PDB_BEG])-1,
                'aend' => intval($region[XML::ATTR_PDB_END])-1,
            ];
        }
        return $ret;
    }
}
