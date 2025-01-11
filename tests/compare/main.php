<?php

require "vendor/autoload.php";

use Tmdet\Compare\Xml\Xml3;
use Tmdet\Compare\Xml\Xml4;
use Tmdet\Compare\Compare;

const PDBTM_3_BASEPATH = "/home/tusi/works/pdbtm_3.0/data/database";
const PDBTM_4_BASEPATH = "/bigdisk/users/tusi/pdbtm/data/database";

function filePath(string $id, string $basePath) : string {
    return $basePath."/".$id[1].$id[2]."/".$id.".xml";
}

foreach(file($argv[1]) as $line) {
    list($id) = explode(' ',trim($line));
    $xml3 = new Xml3(filePath($id,PDBTM_3_BASEPATH));
    $xml4 = new Xml4(filePath($id,PDBTM_4_BASEPATH));
    $compare = new Compare($xml3,$xml4);
    $tmpStatus = $compare->tmpStatus();
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
        echo $id." FAILED ".$tmpStatus[1]." - ".$tmpStatus[2]."\n";
    }
}
