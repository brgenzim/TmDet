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

    public static array $allowedNonStandardResidues = [
        //"ACE", "SO4", "HOH"
        '0G6', '0HH', '12V', '1PE', '1SY', '1YD', '2BA', '3PH', '3WL', '5GP', '6R2', '9NO', '9NX', '9O9', '9OC', '9RL', 'ACP', 'ACT', 'ADN', 'ADP',
        'AJI', 'AKG', 'ALA', 'AMP', 'ANP', 'APA', 'ARG', 'ASN', 'ASP', 'ATP', 'AZA', 'AZI', 'BCT', 'BDP', 'BEN', 'BGC', 'BIS', 'BMA', 'BME',  'BR',
        'BSP', 'BT6', 'C2E',  'CA',  'CD', 'CDL', 'CHT', 'CK9',  'CL', 'CNZ',  'CO', 'CO3', 'CO8', 'CSO', 'CYS', 'D7N', 'DIH', 'DMS', 'DOD', 'DTT',
        'DTV', 'E64', 'EAA', 'ECH', 'EDO', 'EMC', 'EMT', 'EPE', 'FAD',  'FE', 'FE2', 'FES', 'FLC', 'FMN', 'FOM', 'FUA', 'FUC', 'GAL', 'GCP', 'GDN',
        'GDP', 'GEP', 'GLC', 'GLN', 'GLU', 'GLY', 'GNP', 'GOL', 'GSH', 'GTS', 'GTX', 'H34', 'HEC', 'HED', 'HEM',  'HG', 'HIS', 'HOH', 'ILE', 'IMD',
        'IOD', 'IP9', 'IPF', 'JJS', 'JJV', 'JZR',   'K', 'KLO', 'KML', 'KMU', 'KNM', 'KOX', 'KPO', 'KS9', 'KT6', 'KUL', 'KUS', 'LEU',  'LI', 'LYS',
        'M1A', 'MAN', 'MES', 'MET',  'MG', 'MLC', 'MLT', 'MMA',  'MN', 'MPD', 'MRD', 'MYR',  'NA', 'NAG', 'NAI', 'NAP', 'NCA', 'NDG', 'NDP', 'NH4',
         'NI', 'NOX', 'NPO', 'OXM', 'PCT', 'PCW', 'PEE', 'PEG', 'PG4', 'PGE', 'PHE', 'PLM', 'PLP', 'PMP', 'PO4', 'POL', 'PP3', 'PPK', 'PRO', 'PTQ',
        'PXY', 'QTH', 'QTK', 'QTT', 'QU5', 'QWH', 'QWN', 'QXW', 'ROQ', 'SAH', 'SAM', 'SDQ', 'SER', 'SF4', 'SMA', 'SO4', 'STR', 'TAM', 'THR', 'TMO',
        'TRP', 'TRS', 'TYR', 'U2G', 'UNL', 'UNU', 'UNX', 'UPA',  'UQ', 'VAL', 'XL3', 'XPE',  'ZN',
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
        Assert::assertFileExists($this->entFile);
        Assert::assertFileExists($this->cifFile);

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
                if (!in_array($residue, static::$standardResidues)
                    && !in_array($residue, static::$allowedNonStandardResidues)) {

                    throw new InvalidArgumentException(
                        "{$this->pdbCode} ($chainName) has non-standard residue: '$residue'"
                    );
                }
            }
        }
    }
}
