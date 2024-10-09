<?php

namespace Unitmp\TmdetTest\Utils;

use Unitmp\TmdetTest\Exceptions\FileException;

class PdbEnt {
    protected array $polymerChains = [];
    protected array $compoundLines = [];
    protected array $seqResLines = [];
    protected bool $invalidChainList = true;
    const CHAIN_COLUMN = 12;

    public function parse(string $fileName): static {

        $this->invalidChainList = true;

        $pdbLines = explode(PHP_EOL, static::readGzipFile($fileName));
        $this->compoundLines = preg_grep('/^COMPND.+(MOL_ID|MOLECULE|CHAIN):.+/', $pdbLines);
        $this->seqResLines = preg_grep('/^SEQRES.+/', $pdbLines);

        return $this;
    }

    public function getPloymerChains(): array|false {
        if (!$this->invalidChainList) {
            return $this->polymerChains;
        }

        $polymerChainsSeqRes = [];
        foreach ($this->seqResLines as $line) {
            $chain = $line[static::CHAIN_COLUMN - 1];
            // no checking logic, just simple assignment
            $polymerChainsSeqRes[$chain] = $chain;
        }

        if (empty($polymerChainsSeqRes)) {
            return false;
        }

        $moleculeName = '';
        $resultChains = [];
        foreach ($this->compoundLines as $line) {
            $columns = preg_split(
                pattern: '/[ ,;:]/',
                subject: $line,
                flags: PREG_SPLIT_NO_EMPTY
            );
            // Process only the first molecule
            if ($columns[1] == 'MOL_ID' && $columns[2] != '1') {
                break;
            }
            // The lines are in pairs: MOLECULE and CHAIN
            // Chains follows molecule information,
            // but this check is applied before molecule check.
            // (simplicity)
            if ($columns[2] == 'CHAIN') {
                // the first char is the chain name
                $chains = array_slice($columns, 3);
                if (!empty(array_intersect($chains, $polymerChainsSeqRes))) {
                    $chains = array_fill_keys($chains, $moleculeName);
                    $resultChains = array_merge($resultChains, $chains);
                }
                continue;
            }

            $matches = [];
            if (preg_match('/^.+MOLECULE: (.+);/', $line, $matches)) {
                $moleculeName = $matches[1];
            }
        }

        $this->polymerChains = $resultChains;
        $this->invalidChainList = false;
        return $this->polymerChains;
    }


    public static function readGzipFile($filePath): string {
        if (!file_exists($filePath) || !is_readable($filePath)) {
            throw new FileException("File not found or not readable: $filePath");
        }

        $file = gzopen($filePath, 'rb');
        if (!$file) {
            throw new FileException("Failed to open gzip file: $filePath");
        }

        $content = '';
        while (!gzeof($file)) {
            // Read 4KB at a time
            $content .= gzread($file, 4096);
        }

        // Close the file
        gzclose($file);

        return $content;
    }

}