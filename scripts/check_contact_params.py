"""Check RMD contact parameters and boundary_penetration."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
from src.importers.rmd_lexer import RMDLexer

MODELS = os.path.join(os.path.dirname(__file__), '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')

entities = RMDLexer().load(RMD)
for ent in entities:
    if ent.kind == 'CONTACT':
        for k, v in sorted(ent.attrs.items()):
            print(f'  {k}={v}')
