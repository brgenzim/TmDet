<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Services\DsspRunner;
use Unitmp\TmdetTest\Services\PdbtmDsspRunner;
use Unitmp\TmdetTest\Services\TmDetDsspRunner;

class TmDetDsspTest extends TestCase {

    const PDB_ZFS_DIR = '/zfs/databases/UniTmp/PDB/data/structures/divided/mmCIF';
    const PDB_ENT_ZFS_DIR = '/zfs/databases/UniTmp/PDB/data/structures/divided/pdb';

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
    #[DataProvider('getCifPaths')]
    public function test_whole_archive(string $cifPath) {

        // Pre-check
        if (static::isStructureUnsupported($cifPath)) {
            $this->markTestSkipped('DNA/RNA structures are unsupported');
            return;
        }

        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner($cifPath);
        $dsspRunner = DsspRunner::createRunner($cifPath);
        // Act
        $dsspSuccess = $dsspRunner->exec();
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
        $this->assertEquals($dsspRunner->dssps, $tmdetRunner->dssps);
    }

    //
    // Util functions
    //

    /**
     * Data provider for test_whole_archive
     */
    public static function getCifPaths(): array {

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
        return array_slice($files, 4500, 20);
    }

    public static function isStructureUnsupported(string $cifPath): bool {
        exec("zgrep -P 'DNA|RNA' '$cifPath'", $lines, $resultCode);
        return $resultCode === 0;
    }
}
