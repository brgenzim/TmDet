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
foreach(file($argv[1]) as $code) {
    $code = trim($code);
    $xml4 = new Xml4(filePath($code,PDBTM_4_BASEPATH));
    $test = true;
    if ($xml4->getTmpStatus() != "yes") {
        $test = false;
    }
    echo "$code:".($test?"OK":"Not OK")."\n";
}
