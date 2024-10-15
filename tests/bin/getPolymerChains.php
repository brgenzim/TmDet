<?php

use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Utils\PdbEnt;

require_once 'vendor/autoload.php';

$code = $argv[1];

$subDir = '/' . $code[1] . $code[2] . '/';
$fileName = "pdb$code.ent.gz";
$entPath = FileSystem::PDB_ENT_ZFS_DIR . $subDir . $fileName;


if (!file_exists($entPath)) {
    throw new RuntimeException("PDB file of '$code' not found '$entPath'");
}

$polymerChains = (new PdbEnt())->parse($entPath)->getPloymerChains();

echo json_encode($polymerChains, JSON_PRETTY_PRINT|JSON_OBJECT_AS_ARRAY|JSON_UNESCAPED_SLASHES) . PHP_EOL;
