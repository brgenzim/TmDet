<?php

require "vendor/autoload.php";

use Tmdet\Compare\Xml\Xml3;
use Tmdet\Compare\Xml\Xml4;
use Tmdet\Compare\Compare;

const AF_ROOT = "/bigdisk/databases/AlphaFold/data";
const TMDET_ROOT = "/bigdisk/users/tusi/pdbtm/data/AlphaFold";

function filePath(string $id, string $basePath) : string {
    return $basePath."/".$id[1].$id[2]."/".$id.".xml";
}

if ($argc != 2) {
    return "Usage: php Main.php ProteomId";
}
$idsFile = AF_ROOT."/".$argv[1]."/ids.txt";
if (!file_exists($idsFile)) {
    return "Could not find ids file: $idsFile.";
}
foreach(file($idsFile) as $line) {
    $id = trim($line);
    $xml3 = new Xml3(TMDET_ROOT."/old/".$argv[1]."/".$id."_tmdet.xml");
    $xml4 = new Xml4(TMDET_ROOT."/new/".$argv[1]."/".$id.".xml");
    $compare = new Compare($xml3,$xml4);
    $tmpStatus = $compare->tmpStatus();
    if ($tmpStatus[1] == "missing") {
        if ($tmpStatus[1] = "no");
    }
    $tmpStatus[0] = ($tmpStatus[1] == $tmpStatus[2]);
    if ($tmpStatus[0])  {
        if ($tmpStatus[1] == "no") {
            echo $id."  no OK\n";
        }
        else {
            $proteinType = $compare->proteinType();
            if ($proteinType[0]) {
                $chains = $compare->chains();
                if (empty($chains)) {
                    echo $id." yes OK\n";
                }
                else {
                    echo $id." FAILED\n\t".implode("\n\t",$chains)."\n";
                }
            }
            else {
                echo $id." FAILED ".$proteinType[1]." - ".$proteinType[2]."\n";
            }
        }
    }
    else {
        if ( $tmpStatus[1] == "missing") {
            echo $id."  no OK\n";
        }
        else {
            echo $id." FAILED ".$tmpStatus[1]." - ".$tmpStatus[2]."\n";
        }
    }
}
