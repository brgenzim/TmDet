<?php

namespace Unitmp\TmdetTest\Services;

use RuntimeException;
use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Utils\PdbEnt;

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

    public array $polymerChains = [];
    public array $dssps = [];
    public string $newTmdetFile = '';
    public string $oldTmdetFile = '';
    public array $oldData = [];
    public array $newData = [];
    public string $entFile = '';
    public string $pdbCode = '';
    public string $tmdetLogFile = '';
    public bool $enableOverwrite = false;
    // Enables tmdet -no_forcedel option.
    public bool $noForceDelete = false;

    public function __construct(string $execPath, array $commandParams, string $entFile) {
        parent::__construct($execPath, $commandParams);
        $this->disableOutputParsing = true;
        $this->entFile = $entFile;
        $this->polymerChains = (new PdbEnt())->parse($entFile)->getPloymerChains();
        $this->pdbCode = parent::getPdbCodeFromPath($entFile);
        if (!file_exists(static::TMDET_DB_PATH)) {
            mkdir(static::TMDET_DB_PATH);
        }
        // output artifacts
        $dbDir = static::TMDET_DB_PATH;
        $subDir = substr($this->pdbCode, 1, 2);
        $this->newTmdetFile = "$dbDir/$subDir/$this->pdbCode.xml";
        $this->tmdetLogFile = "$dbDir/$subDir/$this->pdbCode.tmdet.log";
        $this->oldTmdetFile = FileSystem::PDBTM30_DATA_DIR . "/database/$subDir/{$this->pdbCode}.xml";

        // $oldData = static::readPdbtmXml($this->oldTmdetFile);
        // $this->oldData = $this->processRegions($oldData);
        $this->oldData = static::readPdbtmXml($this->oldTmdetFile);
    }

    public function exec(): bool {

        if (!file_exists($this->oldTmdetFile)) {
            throw new RuntimeException("Old Tmdet XML of {$this->pdbCode} not found '{$this->oldTmdetFile}'");
        }

        if (file_exists($this->newTmdetFile) && !$this->enableOverwrite) {
            printf("$this->pdbCode: tmdet exec skipped, since output file exists: "
                . "{$this->newTmdetFile} - enableOverwrite suppresses this check\n");
            $result = true;
        } else {
            $result = parent::exec();
            // We do not filter and parse in this kind of runner,
            // but save the output into logfile.
            $output = implode(PHP_EOL, $this->outputLines);
            file_put_contents($this->tmdetLogFile, $output);

            if (!file_exists($this->newTmdetFile)) {
                throw new RuntimeException("New Tmdet XML of {$this->pdbCode} not found '{$this->newTmdetFile}'");
            }
        }

        // TODO: remove it later - if processRegions is not needed
        // $newData = static::readPdbtmXml($this->newTmdetFile);
        // $this->newData = $this->processRegions($newData);

        $this->newData = static::readPdbtmXml($this->newTmdetFile);

        return $result;
    }

    // protected function processRegions(array $pdbtmContent): array {
    //     $regions = [];
    //     foreach ($pdbtmContent['chains'] as $chainId => $chainItem) {
    //         $typeSequence = '';
    //         foreach ($chainItem['regions'] as $region) {
    //             $count = $region['seq_end'] - $region['seq_beg'] + 1;
    //             $typeSequence .= str_repeat($region['type'], $count);
    //             // SHORT region strings: $typeSequence .= $region['type'];
    //         }
    //         $regions[$chainId] = $typeSequence;
    //     }
    //     return [
    //         'tmType' => $pdbtmContent['tmType'],
    //         'isTmp' => $pdbtmContent['isTmp'],
    //         'chains' => $regions
    //     ];
    // }

    protected function filterOutputLines(array $lines): array {
        // output will not be parsed
        return [];
    }

    protected function parseOutputLine(string $line): string {
        // no operation
        return $line;
    }

    public static function createRunner(string $entFile, bool $noForceDelete = false): static {

        $pdbCode = static::getPdbCodeFromPath($entFile);
        $params = [
            "-p=$pdbCode",
            '-create'
        ];
        if ($noForceDelete) {
            $params[] = '-no_forcedel';
        }
        $runner = new Tmdet30Runner(static::EXEC, $params, $entFile);
        $runner->noForceDelete = $noForceDelete;
        return $runner;
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
            $num_tm = (int) $chain['NUM_TM'];
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
        $deletedChains = [];
        $addedChains = [];
        if (isset($xml->BIOMATRIX)) {
            foreach ($xml->BIOMATRIX->DELETE as $deleteItem) {
                $deletedChains[] = (string) $deleteItem['CHAINID'];
            }
            sort($deletedChains);

            foreach ($xml->BIOMATRIX->MATRIX as $matrix) {
                foreach ($matrix->APPLY_TO_CHAIN as $applyItem) {
                    $addedChains[] = (string) $applyItem['NEW_CHAINID'];
                }
            }
            sort($addedChains);
        }

        return [
            'code' => (string)$xml['ID'],
            'chains' => $chains,
            'deletedChains' => $deletedChains,
            'addedChains' => $addedChains,
            'isTmp' => ((string)$xml['TMP'] === 'yes') ? true : false,
            'tmType' => (string)$xml->RAWRES->TMTYPE
        ];
    }

    public static function getPdbtmCodes(string $xmlOfAllEntries): array {
        $root = simplexml_load_file($xmlOfAllEntries);

        $pdbtmProteins = [];
        foreach ($root->pdbtm as $item) {
            $code = (string)$item['ID'];
            $pdbtmProteins[$code] = [
                'code' => $code,
                'isTmp' => ((string)$item['TMP'] === 'yes') ? true : false,
                'tmType' => (string)$item->RAWRES->TMTYPE
            ];
        };
        return $pdbtmProteins;
    }

    public function getCommandLine(): string {
        return $this->commandLine;
    }
}
