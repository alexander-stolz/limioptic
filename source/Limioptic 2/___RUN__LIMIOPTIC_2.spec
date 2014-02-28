# -*- mode: python -*-
a = Analysis(['___RUN__LIMIOPTIC_2.py'],
             pathex=['G:\\GitHub\\limioptic\\source\\Limioptic 2'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='___RUN__LIMIOPTIC_2.exe',
          debug=False,
          strip=None,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='___RUN__LIMIOPTIC_2')
