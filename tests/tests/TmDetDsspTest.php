<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\DsspRunner;
use Unitmp\TmdetTest\Services\PdbtmDsspRunner;
use Unitmp\TmdetTest\Services\StructurePreFilter;
use Unitmp\TmdetTest\Services\TmDetDsspRunner;

class TmDetDsspTest extends TestCase {

    const PDB_ZFS_DIR = FileSystem::PDB_ZFS_DIR;
    const PDB_ENT_ZFS_DIR = FileSystem::PDB_ENT_ZFS_DIR;

    #[Group('dssp')]
    public function test_1afo_dssp() {
        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner('1AFO');
        $dsspRunner = DsspRunner::createRunner('1AFO');
        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        $dsspSuccess = $dsspRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
    }

    #[Group('dssp')]
    public function test_pdbtm_dssp_integrity_with_132l() {
        // Arrange
        $runner = PdbtmDsspRunner::createRunner(static::PDB_ZFS_DIR . '/32/pdb132l.ent.gz');
        // Act
        $isSuccess = $runner->exec();
        // Assert
        $this->assertTrue($isSuccess);
        $this->assertEquals('A', $runner->chains[0]);
    }

    #[Group('dssp')]
    public function test_pdbtm_dssp_integrity_with_1aft() {
        $path = static::PDB_ZFS_DIR . '/af/1aft.cif.gz';
        $filter = new StructurePreFilter($path);
        $this->assertEquals('1aft', $filter->pdbCode);
        $filter->checkEntryFiles();
    }

    #[Group('dssp')]
    public function test_pdbtm_dssp_integrity_with_7rh7_chain_w() {
        // Arrange
        $pdbtmRunner = PdbtmDsspRunner::createRunner(static::PDB_ENT_ZFS_DIR . '/rh/pdb7rh7.ent.gz');
        $tmdetRunner = TmDetDsspRunner::createRunner(static::PDB_ZFS_DIR . '/rh/7rh7.cif.gz');
        // Act
        $isSuccess = $pdbtmRunner->exec();
        $isTmDetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($isSuccess);
        $this->assertTrue($isTmDetSuccess);
        $this->assertEquals($pdbtmRunner->dssps['W'], $tmdetRunner->dssps['W']);
        $this->assertEquals($pdbtmRunner->dssps, $tmdetRunner->dssps);
    }

    #[Group('dssp')]
    public function test_7rh7_dssp() {
        // Arrange
        $pdbtmRunner = PdbtmDsspRunner::createRunner(static::PDB_ENT_ZFS_DIR . '/rh/pdb7rh7.ent.gz');
        $tmdetRunner = TmDetDsspRunner::createRunner(static::PDB_ZFS_DIR . '/rh/7rh7.cif.gz');
        // Act
        $isSuccess = $pdbtmRunner->exec();
        $isTmDetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($isSuccess);
        $this->assertTrue($isTmDetSuccess);
        $this->assertEquals($pdbtmRunner->dssps, $tmdetRunner->dssps);
    }

    #[Group('dssp')]
    public function test_8a1e_dssp() {
        // Arrange
        $pdbtmRunner = PdbtmDsspRunner::createRunner(static::PDB_ENT_ZFS_DIR . '/a1/pdb8a1e.ent.gz');
        $tmdetRunner = TmDetDsspRunner::createRunner(static::PDB_ZFS_DIR . '/a1/8a1e.cif.gz');
        // Act
        $isSuccess = $pdbtmRunner->exec();
        $isTmDetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($isSuccess);
        $this->assertTrue($isTmDetSuccess);
        $this->assertEquals($pdbtmRunner->dssps, $tmdetRunner->dssps);
    }

    #[Group('dssp')]
    public function test_6a1t_dssp() {
        // Arrange
        $pdbtmRunner = PdbtmDsspRunner::createRunner(static::PDB_ENT_ZFS_DIR . '/a1/pdb6a1t.ent.gz');
        $tmdetRunner = TmDetDsspRunner::createRunner(static::PDB_ZFS_DIR . '/a1/6a1t.cif.gz');
        // Act
        $isSuccess = $pdbtmRunner->exec();
        $isTmDetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($isSuccess);
        $this->assertTrue($isTmDetSuccess);
        $this->assertEquals($pdbtmRunner->dssps, $tmdetRunner->dssps);
    }

    #[Group('dssp')]
    #[DataProvider('getCifPaths')]
    //#[DataProvider('staticCifPathProvider')]
    public function test_whole_archive(string $cifPath) {

        // Pre-check
        $filter = new StructurePreFilter($cifPath);
        $exception = null;
        if (!$filter->tryCheckEntryFiles($exception)) {
            printf("Skipping test: %s\n", $exception->getMessage());
            $this->markTestSkipped('Unsupported CIF-ENT pair ' . $cifPath);
        }

        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner($cifPath);
        $dsspRunner = PdbtmDsspRunner::createRunner($filter->entFile);
        // Act
        $dsspSuccess = $dsspRunner->exec();
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        // $this->assertTrue($tmdetSuccess);
        // $this->assertTrue($dsspSuccess);
        $this->assertEquals($dsspRunner->dssps, $tmdetRunner->dssps);
    }

    //
    // Util functions
    //

    /**
     * Data provider for test_whole_archive
     */
    public static function getCifPaths(): array {

        printf("Scanning PDB directories...");
        $dirs = scandir(static::PDB_ZFS_DIR);
        if (!$dirs || count($dirs) == 2) {
            return [];
        }

        $files = [];
        foreach (array_slice($dirs, 2) as $directory) {
            $directory = static::PDB_ZFS_DIR . "/$directory";
            $items = scandir($directory);
            $items = preg_grep('/^.+cif.gz$/', $items);
            foreach($items as $item) {
                $files[ $item ] = [ "$directory/$item" ];
            }
        }
        // TODO: remove slice later
        //return array_slice($files, 100, 25);
        //return array_slice($files, 100, 1000);
        return array_slice($files, 40000, 1000);
        //return $files;
    }

    public static function staticCifPathProvider(): array {
        return [
            '154l.cif.gz' => [ static::PDB_ZFS_DIR . '/54/154l.cif.gz' ],
            '1a0h.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/1a0h.cif.gz' ],
            '1a0t.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/1a0t.cif.gz' ],
            '2a0z.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/2a0z.cif.gz' ],
            '3a09.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/3a09.cif.gz' ],
            '3a0e.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/3a0e.cif.gz' ],
            '4a05.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/4a05.cif.gz' ],
            '4a0p.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/4a0p.cif.gz' ],
            '5a05.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/5a05.cif.gz' ],
            '5a09.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/5a09.cif.gz' ],
            '5a0a.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/5a0a.cif.gz' ],
            '5a0b.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/5a0b.cif.gz' ],
            '5a0c.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/5a0c.cif.gz' ],
            '6a0j.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/6a0j.cif.gz' ],
            '6a0k.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/6a0k.cif.gz' ],
            '6a0l.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/6a0l.cif.gz' ],
            '7a0k.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/7a0k.cif.gz' ],
            '7a0q.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/7a0q.cif.gz' ],
            '7a0t.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/7a0t.cif.gz' ],
            '8a0y.cif.gz' => [ static::PDB_ZFS_DIR . '/a0/8a0y.cif.gz' ],
            '1a14.cif.gz' => [ static::PDB_ZFS_DIR . '/a1/1a14.cif.gz' ],
            '6a1t.cif.gz' => [ static::PDB_ZFS_DIR . '/a1/6a1t.cif.gz' ],
            '8a16.cif.gz' => [ static::PDB_ZFS_DIR . '/a1/8a16.cif.gz' ],
            '8a17.cif.gz' => [ static::PDB_ZFS_DIR . '/a1/8a17.cif.gz' ],
            '8a1f.cif.gz' => [ static::PDB_ZFS_DIR . '/a1/8a1f.cif.gz' ],
        ];
    }

    public static function isStructureUnsupported(string $cifPath): bool {
        exec("zgrep -P 'DNA|RNA' '$cifPath'", $lines, $resultCode);
        return $resultCode === 0;
    }
}
