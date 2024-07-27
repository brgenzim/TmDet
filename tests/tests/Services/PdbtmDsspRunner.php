<?php

namespace Unitmp\TmdetTest\Services;

use RuntimeException;

class PdbtmDsspRunner extends AbstractProcessRunner {

    const PDB_CACHE_PATH = '/tmp';
    const LD_LIBRARY_PATH = '/home/A/csongor/dev/PdbLib.C__stale/lib';
    const EXEC = 'export PDB_CACHE_PATH=' . PdbtmDsspRunner::PDB_CACHE_PATH
        . '; export LD_LIBRARY_PATH=' . PdbtmDsspRunner::LD_LIBRARY_PATH
        // . '; /home/tusi/works/pdbtm_3.0/TmdetUtils/bin/dssp';
        . '; /home/A/csongor/dev/TmdetUtils/bin/dssp';
    const DUMP_HEADER_LINE_PREFIX = 'ShowProt:';
    const CHAIN_COLUMN = 1;
    const STRUCTURE_COLUMN = 5;

    public array $chains = [];
    public array $dssps = [];
    public string $dsspCacheFile = '';
    public string $entFile = '';
    public string $pdbCode = '';

    public function __construct(string $execPath, array $commandParams, string $entFile) {
        parent::__construct($execPath, $commandParams);
        $this->entFile = $entFile;
        $this->pdbCode = parent::getPdbCodeFromPath($entFile);
    }

    protected function filterOutputLines(array $lines): array {

        $cacheDir = static::PDB_CACHE_PATH;
        $dsspFiles = explode("\n", trim(`ls -t $cacheDir/$this->pdbCode.*.dssp`));
        var_dump($dsspFiles);
        if (!empty($dsspFiles) && file_exists($dsspFile = array_shift($dsspFiles))) {
            $this->dsspCacheFile = $dsspFile;

            return explode("\n", file_get_contents($this->dsspCacheFile));
        }

        throw new RuntimeException("DSSP cache file of {$this->pdbCode} not found");
    }

    protected function parseOutputLine(string $line): string {

        if (!empty($this->dsspCacheFile)) {
            return $this->parseDsspLine($line);
        }

        return $line;
    }

    protected function parseDsspLine(string $line): string {
        // parse a chain line
        if (preg_match('/^(\w+)_(\w+)/', $line, $matches)) {
            $chain = $matches[2];
            if (!in_array($chain, $this->chains)) {
                $this->chains[] = $chain;
            }
        } elseif (count(($columns = explode(' ', $line))) == 7) {
            $column = $columns[static::CHAIN_COLUMN - 1];
            if (!array_key_exists($column, $this->dssps)) {
                $this->dssps[$column] = '';
            }
            $secondaryStructure = $columns[static::STRUCTURE_COLUMN - 1];
            if ($secondaryStructure == ' ') {
                $secondaryStructure = '-';
            }
            $this->dssps[$column] = $this->dssps[$column]
                . $secondaryStructure;
        }
        return $line;
    }

    public static function createRunner(string $entFile): static {

        $params = [
            '-i', "'$entFile'",
            '--output-format', 'dssp'
        ];
        return new PdbtmDsspRunner(static::EXEC, $params, $entFile);
    }
}
