<?php

namespace Tmdet\Compare;

use Tmdet\Compare\Xml\Xml3;
use Tmdet\Compare\Xml\Xml4;

class Compare {
    private $xml3;
    private $xml4;

    public function __construct(Xml3 $xml3, Xml4 $xml4) {
        $this->xml3 = $xml3;
        $this->xml4 = $xml4;
    }

    public function tmpStatus() : array {
        $ret3 = $this->xml3->getTmpStatus();
        $ret4 = $this->xml4->getTmpStatus();
        return  [ $ret3 == $ret4, $ret3, $ret4];
    }

    public function proteinType() : array {
        $ret3 = $this->xml3->getProteinType();
        $ret3 = ($ret3 == "Ca_Tm"?"Tm_Alpha":$ret3);
        $ret4 = $this->xml4->getProteinType();
        $ret4 = ($ret4=="Tm_Mixed"?$ret3:$ret4);
        return [ $ret3 == $ret4, $ret3, $ret4];
    }

    public function chains() : array {
        $ret = [];
        $ret3 = $this->xml3->getChains();
        $ret4 = $this->xml4->getChains();
        foreach($ret3 as $id => $chainData) {
            if (!isset($ret4[$id])) {
                //$ret[] = "Missing chain: $id";
            }
            else {
                $ret = array_merge($ret,$this->chain($ret3[$id],$ret4[$id]));
            }
        }
        return $ret;
    }

    private function chain(array $ch3, array $ch4) : array {
        $ret = [];
        if ( abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 0 
            && abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 1
            && abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 2
            && abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 4) {
            $ret[] = "Number of tms is different for chain ".$ch3['id'].": ".$ch3['numtm']." vs ".$ch4['numtm'];
        }
        /*if ($ch3['seq'] != $ch4['seq']) {
            $ret[] = "Sequence is different for chain ".$ch3['id'].": ".$ch3['seq']." vs ".$ch4['seq'];
        }*/
        if (empty($ret)) {
            $ret = $this->regions($ch3['regions'],$ch4['regions']);
        }
        return $ret;
    }

    private function regions(array $r3, array $r4) : array {
        $ret = [];

        return $ret;
    }
}
