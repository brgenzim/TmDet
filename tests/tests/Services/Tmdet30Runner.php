<?php

namespace Unitmp\TmdetTest\Services;

use RuntimeException;
use Unitmp\TmdetTest\Constants\FileSystem;

class Tmdet30Runner extends AbstractProcessRunner {

    // export TMDET_DATABASE_PATH="$TARGET_DIR"
    // export PDBINFO_LINES="HEADER:TITLE :COMPND:SOURCE:SITE "
    // export PDBMODE="subdirs"
    // export PDBPATH="/data/wwPDB/data/structures/divided/pdb/"
    // export PDBPATH="$TMDET_DATABASE_PATH"
    // export PDBMODELPATH="$TMDET_DATABASE_PATH"
    // export PDBTYPE="compressed"
    // export PDB_ENT="ent"
    // export PDB_INF="HEADER:TITLE :COMPND:SOURCE:SITE:REMARK"

    // export PDB_CACHE_PATH="$TMDET_DATABASE_PATH"



    const PDB_CACHE_PATH = '/tmp/pdb-cache-tmdet-3.0';
    const EXEC = 'export PDB_CACHE_PATH=' . Tmdet30Runner::PDB_CACHE_PATH
        . '; export PDBMODE="subdirs"'
        . '; export TMDET_DATABASE_PATH=' . Tmdet30Runner::PDB_CACHE_PATH
        . '; export PDBPATH=' . FileSystem::PDBTM30_DATA_DIR . '/database'
        . '; /home/tusi/works/pdbtm_3.0/TmdetUtils/bin/tmdet';

    public array $chains = [];
    public array $dssps = [];
    public string $tmdetFile = '';
    public string $entFile = '';
    public string $pdbCode = '';

    public function __construct(string $execPath, array $commandParams, string $entFile) {
        parent::__construct($execPath, $commandParams);
        $this->entFile = $entFile;
        $this->pdbCode = parent::getPdbCodeFromPath($entFile);
        if (!file_exists(static::PDB_CACHE_PATH)) {
            mkdir(static::PDB_CACHE_PATH);
        }
    }

    protected function filterOutputLines(array $lines): array {

        $cacheDir = static::PDB_CACHE_PATH;
        $subDir = substr($this->pdbCode, 1, 2);
        $tmdetFile = "$cacheDir/$subDir/$this->pdbCode.xml";
        if (file_exists($tmdetFile)) {
            $this->tmdetFile = $tmdetFile;

            return explode("\n", file_get_contents($this->tmdetFile));
        }

        throw new RuntimeException("Tmdet XML of {$this->pdbCode} not found '{$this->tmdetFile}'");
    }

    protected function parseOutputLine(string $line): string {
        return $line;
    }

    public static function createRunner(string $entFile): static {

        $pdbCode = static::getPdbCodeFromPath($entFile);
        $params = [
            "-p=$pdbCode",
            '-create'
        ];
        return new Tmdet30Runner(static::EXEC, $params, $entFile);
    }
}
