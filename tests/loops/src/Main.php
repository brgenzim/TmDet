<?php

require "vendor/autoload.php";

use Tmdet\Compare\Xml\Xml4;

const PDBTM_4_BASEPATH = "/bigdisk/users/tusi/pdbtm/data/database";

function filePath(string $id, string $basePath) : string {
    return $basePath."/".$id[1].$id[2]."/".$id.".xml";
}

if ($argc != 2) {
    return "Usage: php Main.php data/codes.txt";
}
if (!file_exists($argv[1])) {
    return "Input file does not exists.";
}
foreach(file($argv[1]) as $line) {
    
    list($code) = explode(' ',trim($line));
    $xml4 = new Xml4(filePath($code,PDBTM_4_BASEPATH));
    if ($xml4->getTmpStatus()) {
        echo "$code:".($xml4->hasRequestedRegion('L')||$xml4->hasRequestedRegion('K')?"OK":"Not OK")."\n";
    }
    else {
        echo "$code: ---NOT TMP---\n";
    }
}
