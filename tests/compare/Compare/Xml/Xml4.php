<?php

namespace Tmdet\Compare\Xml;

use SimpleXMLElement;
use Tmdet\Compare\Helper;
use Tmdet\Compare\Xml\Constant4 as XML;

class Xml4 {
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
        return $this->xml ? $this->xml[XML::ATTR_TRANSMEMBRANE]->__toString() : "missing";
    }

    public function getProteinType() {
        return $this->xml ? $this->xml->{XML::NODE_RAWDATA}->{ Xml::NODE_TMTYPE }->__toString() : "missing";
    }

    public function getChains() : array {
        $ret = [];
        if ($this->xml !== false) {
            foreach($this->xml->{XML::NODE_CHAINS}->{XML::NODE_CHAIN} as $chain) {
                $id = $chain[XML::ATTR_AUTH_ID]->__toString();
                $id = ($id=="_"?"A":$id);
                $type = $chain[XML::ATTR_TYPE]->__toString();
                if ($type == 'alpha' || $type == 'beta') {
                    $ret[$id] = [
                        'id' => $id,
                        'type' => $type,
                        'numtm' => $chain[XML::ATTR_NUM_TM]->__toString(),
                        'seq' => Helper::clearSeq($chain->{XML::NODE_SEQENCE}->__toString()),
                        'regions' => $this->getRegions($chain),
                    ];
                }
            }
        }
        return $ret;
    }

    private function getRegions(SimpleXMLElement $chainNode) : array {
        $ret = [];
        foreach($chainNode->{XML::NODE_REGIONS}->{XML::NODE_REGION} as $region) {
            $ret[] = [
                'start' => $region[XML::ATTR_START_LABEL_ID],
                'end' => $region[XML::ATTR_END_LABEL_ID],
                'type' => $region[XML::ATTR_TYPE]
            ];
        }
        return $ret;
    }
}
