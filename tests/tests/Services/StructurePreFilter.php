<?php

namespace Unitmp\TmdetTest\Services;

use InvalidArgumentException;
use PHPUnit\Framework\Assert;
use RuntimeException;
use Unitmp\TmdetTest\Constants\FileSystem;

class StructurePreFilter {

    const GEMMI_CMD = "gemmi residues -m '/1/*' '%s'";
    public readonly array $commandParams;
    public array $outputLines = [];
    public string $pdbCode;
    public readonly string $entFile;
    public readonly string $cifFile;
    public array $entChains;
    public array $cifChains;

    public static array $standardResidues = [
        "ALA", "CYS", "ASP", "GLU",
        "PHE", "GLY", "HIS", "ILE",
        "LYS", "LEU", "MET", "ASN",
        "PRO", "GLN", "ARG", "SER",
        "THR", "VAL", "TRP", "TYR",
    ];

    public function __construct(string $cifPath) {
        $this->pdbCode = AbstractProcessRunner::getPdbCodeFromPath($cifPath);
        $this->cifFile = $cifPath;
        $this->entFile = FileSystem::PDB_ENT_ZFS_DIR
            . '/' . substr($this->pdbCode, 1, 2) // two middle chars of the code
            . "/pdb$this->pdbCode.ent.gz";
    }

    /**
     * Verifies each residue are std. residue in the CIF and ENT files.
     * Compares chains and their residues are identical.
     */
    public function checkEntryFiles(): void {
        $this->cifChains = $this->checkEntryFile($this->cifFile);
        $this->entChains = $this->checkEntryFile($this->entFile);
        Assert::assertEquals($this->entChains, $this->cifChains);
    }

    /**
     * Check residues are standard residues.
     * Returns array of chains with array of residues.
     */
    public function checkEntryFile(string $entryFile): array {
        $lines = null;
        $resultCode = null;

        $command = sprintf(static::GEMMI_CMD, $entryFile);
        exec($command, $lines, $resultCode);

        if ($resultCode !== 0) {
            fprintf(STDERR, "%s: failed: %s\n", $this->pdbCode, $command);
            fprintf(STDERR, "%s: output:\n%s\n", $this->pdbCode, implode("\n", $lines));
            throw new RuntimeException('Failed command: ' . $command);
        }

        $chains = $this->getChains($lines);
        $this->checkChains($chains);

        return $chains;
    }

    protected function getChains(array $lines): array {
        $chains = [];
        foreach ($lines as $line) {
            if (!preg_match('/^(\S+)\s+\d+\s+(\S+)\s+.+/', $line, $matches)) {
                continue;
            }
            $chain = $matches[1];
            $residue = $matches[2];

            if (!in_array($chain, array_keys($chains))) {
                $chains[$chain] = [];
            }

            $chains[$chain][] = $residue;
        }
        return $chains;
    }

    protected function checkChains(array $chains): void {
        foreach ($chains as $chainName => $chain) {
            foreach ($chain as $residue) {
                if (!in_array($residue, static::$standardResidues)) {
                    throw new InvalidArgumentException("{$this->pdbCode} ($chainName) has non-standard residue: '$residue'");
                }
            }
        }
    }
}
