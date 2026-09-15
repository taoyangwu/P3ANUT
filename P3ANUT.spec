# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['P3ANUT.py'],
    pathex=['src'],
    binaries=[],
    datas=[
        ('config.yaml', '.'),
    ],
    hiddenimports=[
        'p3anut_ui',
        'sequenceCounter',
        'multiprocessedPairAssembler',
        'runUnifier',
        'volcanoPlot',
        'upsetPlot',
        'rankingPlot',
        'CLI_VolcanoPlot',
        'CLI_upsetplot',
        'CLI_rankingPlot',
        'utils.FASTA_fileConversion',
        'utils.visualizationGraphs',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='P3ANUT',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
app = BUNDLE(
    exe,
    name='P3ANUT.app',
    icon=None,
    bundle_identifier=None,
)
