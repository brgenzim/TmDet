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
        /*if ( abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 0 
            && abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 1
            && abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 2
            && abs(intval($ch3['numtm']) - intval($ch4['numtm'])) != 4) {
            $ret[] = "Number of tms is different for chain ".$ch3['id'].": ".$ch3['numtm']." vs ".$ch4['numtm'];
        }*/
        if ($ch3['type'] != $ch4['type']) {
            $ret[] = "Type is different for chain ".$ch3['id'].": ".$ch3['type']." vs ".$ch4['type'];
        }
        /*if ($ch3['seq'] != $ch4['seq']) {
            $ret[] = "Sequence is different for chain ".$ch3['id'];
        }*/
        else {
            $res = $this->regions($ch3['regions'],$ch4['regions'],strlen($ch3['seq']));
            if ($res != "") {
                $ret[] = $ch3['id'].": ".$res;
            }
        }
        return $ret;
    }

    private function regions(array $r3, array $r4, int $len) : string {
        $ret = "";
        $fp=$fn=0;
        $top3 = $this->getTopologyAsString($r3,$len);
        $top4 = $this->getTopologyAsString($r4,$len);
        if (strlen($top3) != strlen($top4)) {
            $ret = "Lengths of topology strings are different:\n$top3\n$top4\n";
        }
        elseif ($this->compareTopologies($top3, $top4, $fp, $fn)) {
            $ret= "Topologies are different: fp: $fp fn: $fn\n$top3\n$top4\n";
        }
        return $ret;
    }

    private function getTopologyAsString(array $regions, int $len) : string {
        $top = "";
        for ($i=0; $i<$len; $i++) {
            $top .= '.';
        }
        for($i=0; $i<count($regions); $i++) {
            $reg=$regions[$i];
            if ($reg['type'] == 'H' || $reg['type'] == 'B' || $reg['type'] == 'L') {
                if (($j = intval($reg['start'])) > 0) {
                    $top[$j-1]='.';
                }
                for($j=$reg['start']; $j<=$reg['end']; $j++) {
                    $top[$j] = $reg['type'];
                }
                if ($reg['type'] == 'H'
                    && $i > 0
                    && $i < count($regions)-1
                    && $regions[$i-1]['type'] == $regions[$i+1]['type']
                    && $regions[$i-1]['type'] != 'G'
                    && $regions[$i]['astart']-$regions[$i-1]['aend'] == 1
                    && $regions[$i+1]['astart']-$regions[$i]['aend'] == 1) {
                        $k=(int)(($reg['start']+$reg['end'])/2);
                        $top[$k] = '.';
                }
            }
        }
        
        return $top;
    }

    private function compareTopologies(string $top1, string $top2, int &$fp, int &$fn) : bool {
        $beg=$end=0;
        $ok=$fp=$fn=0;
        $len=strlen($top1);
        while($this->getNextCommonRegion($top1,$top2,$len,$beg,$end)) {
            if ($end-$beg>0) {
                $ok++;
                if ($top1[$beg]!='L') {
                    $this->clearRegion($top1,$beg,$end,$len);
                }
                if ($top2[$beg]!='L') {
                    $this->clearRegion($top2,$beg,$end,$len);
                }
            }
            $beg=$end+1;
        }
        
        $fn=$this->countRegions($top1,$len);
        $fp=$this->countRegions($top2,$len);
        if ($fp>1||$fn>1) {
            echo "$top1<<\n$top2\n";
        }
        return ($fp>1||$fn>1);
    }

    private function getNextCommonRegion(string $top1, string $top2, int $len, int &$beg, int &$end): bool {
        while($beg<$len && ($top1[$beg] == '.' || $top2[$beg] == '.')) {
            $beg++;
        }
        if($beg==$len) {
            return false;
        }
        $end=$beg+1;
        while($end<$len && $top1[$end] != '.' && $top2[$end] != '.') {
            $end++;
        }
        $end--;
        return true;
    }

    private function clearRegion(string &$top, int $beg, int $end, int $len) {
        $i=$beg;
        while($i>=0 && $top[$i] != '.') {
            $top[$i]='.';
            $i--;
        }
        for($i=$beg; $i<=$end; $i++) {
            $top[$i] = '.';
        }
        //if ($i<$len && !$loop) {
            while($i<$len && $top[$i] != '.') {
                $top[$i]='.';
                $i++;
            }
        //}
    }

    private function countRegions(string $top, int $len) {
        $beg=$end=0;
        $ret=0;
        while($this->getNextRegion($top,$len,$beg,$end)) {
            $ret++;
            $beg=$end+1;
        }
        return $ret;
    }

    private function getNextRegion(string $top, int $len, int &$beg, int &$end): bool {
        while($beg<$len && ($top[$beg] == '.' || $top[$beg] == 'L')) {
            $beg++;
        }
        if($beg==$len) {
            return false;
        }
        $end=$beg+1;
        while($end<$len && $top[$end] != '.' && $top[$end] != 'L') {
            $end++;
        }
        $end--;
        return true;
    }
}
