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
        $this->compoundLines = preg_grep('/^COMPND\s+([\d]+)? /', $pdbLines);
        $this->compoundLines = preg_replace('/^COMPND\s+([\d]+)? /', '', $this->compoundLines);
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

        $pattern = '/MOL_ID: (.*?);\s+MOLECULE: (.*?);\s+CHAIN: (.*?);/ms';
        if (empty($polymerChainsSeqRes)
            || preg_match_all($pattern, implode(PHP_EOL, $this->compoundLines), $molecules) == 0) {

            return false;
        }

        $moleculeCount = count($molecules[0]);
        $moleculeName = '';
        $resultChains = [];
        for ($index = 0; $index < $moleculeCount; $index++) {
            // we have 3 groups in the pattern
            // $moleculeId = $molecules[1][$index];
            $moleculeName = $molecules[2][$index];
            $moleculeName = preg_replace("/\s+/s", ' ', $moleculeName);
            $moleculeChains = $molecules[3][$index];
            $moleculeChains = preg_split(
                pattern: '/[ ,;]/',
                subject: $moleculeChains,
                flags: PREG_SPLIT_NO_EMPTY
            );

            // the first char is the chain name
            if (!empty(array_intersect($moleculeChains, $polymerChainsSeqRes))) {
                $chains = array_fill_keys($moleculeChains, $moleculeName);
                $resultChains = array_merge($resultChains, $chains);
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
