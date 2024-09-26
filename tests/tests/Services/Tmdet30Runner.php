<?php

namespace Unitmp\TmdetTest\Services;

use RuntimeException;
use Unitmp\TmdetTest\Constants\FileSystem;

class Tmdet30Runner extends AbstractProcessRunner {

    const TMDET_DB_PATH = '/tmp/tmdet-3.0';
    const EXEC = 'export PDBMODE="subdirs"'
        . '; export TMDET_DATABASE_PATH="' . Tmdet30Runner::TMDET_DB_PATH . '"'
        . '; export PDBPATH="' . FileSystem::PDB_ENT_ZFS_DIR . '"'
        . '; export PDBMODELPATH="$PDBPATH"'
        . '; export PDB_CACHE_PATH="$PDBPATH"'
        . '; export PDB_ENT="ent"'
        . '; export PDBTYPE="compressed"'
        . '; export PDBINFO_LINES="HEADER:TITLE :COMPND:SOURCE:SITE "'
        . '; export PDB_INF="$PDBINFO_LINES"'
        . '; /home/tusi/works/pdbtm_3.0/TmdetUtils/bin/tmdet';

    public array $chains = [];
    public array $dssps = [];
    public string $newTmdetFile = '';
    public string $oldTmdetFile = '';
    public array $oldData = [];
    public array $newData = [];
    public string $entFile = '';
    public string $pdbCode = '';

    public function __construct(string $execPath, array $commandParams, string $entFile) {
        parent::__construct($execPath, $commandParams);
        $this->disableOutputParsing = true;
        $this->entFile = $entFile;
        $this->pdbCode = parent::getPdbCodeFromPath($entFile);
        if (!file_exists(static::TMDET_DB_PATH)) {
            mkdir(static::TMDET_DB_PATH);
        }
        // output artifacts
        $dbDir = static::TMDET_DB_PATH;
        $subDir = substr($this->pdbCode, 1, 2);
        $this->newTmdetFile = "$dbDir/$subDir/$this->pdbCode.xml";
        $this->oldTmdetFile = FileSystem::PDBTM30_DATA_DIR . "/database/$subDir/{$this->pdbCode}.xml";
    }

    public function exec(): bool {

        if (!file_exists($this->oldTmdetFile)) {
            throw new RuntimeException("Old Tmdet XML of {$this->pdbCode} not found '{$this->oldTmdetFile}'");
        }

        $result = parent::exec();

        if (!file_exists($this->newTmdetFile)) {
            throw new RuntimeException("New Tmdet XML of {$this->pdbCode} not found '{$this->newTmdetFile}'");
        }

        $oldData = static::readPdbtmXml($this->oldTmdetFile);
        $this->oldData = $this->processRegions($oldData);

        $newData = static::readPdbtmXml($this->newTmdetFile);
        $this->newData = $this->processRegions($newData);

        return $result;
    }

    protected function processRegions(array $pdbtmContent): array {
        $regions = [];
        foreach ($pdbtmContent as $chainId => $chainItem) {
            $typeSequence = '';
            foreach ($chainItem['regions'] as $region) {
                // $count = $region['seq_end'] - $region['seq_beg'] + 1;
                // $typeSequence .= str_repeat($region['type'], $count);
                $typeSequence .= $region['type'];
            }
            $regions[$chainId] = $typeSequence;
        }
        return $regions;
    }

    protected function filterOutputLines(array $lines): array {
        return [];
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

    // Util function
    public static function readPdbtmXml(string $xmlPath): array {

        $xml = simplexml_load_file($xmlPath);
        if (!$xml) {
            throw new RuntimeException("Pdbtm XML of not found: '$xmlPath'");
        }
        $chains = [];

        foreach ($xml->CHAIN as $chain) {
            $chainID = (string) $chain['CHAINID'];
            $num_tm = (string) $chain['NUM_TM'];
            $type = (string) $chain['TYPE'];
            $regions = [];
            foreach ($chain->REGION as $region) {
                $regions[] = [
                    'seq_beg' => (int) $region['seq_beg'],
                    'pdb_beg' => (int) $region['pdb_beg'],
                    'seq_end' => (int) $region['seq_end'],
                    'pdb_end' => (int) $region['pdb_end'],
                    'type'    => (string) $region['type']
                ];
            }

            $chains[$chainID] = [
                'num_tm' => $num_tm,
                'type' => $type,
                'regions' => $regions
            ];
        }
        return $chains;
    }
}
