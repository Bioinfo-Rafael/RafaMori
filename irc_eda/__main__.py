from __future__ import annotations
import argparse
import json
import sys
from .pipeline import Pipeline,load_config
from .samples import samples


def main(argv=None):
    parser=argparse.ArgumentParser(description='Gaussian IRC → official APOST3D ENPART / easyEDA sobEDA adapters')
    parser.add_argument('command',choices=['inspect','prepare','run-iqa','run-eda','analyze','run-all','check-samples'])
    parser.add_argument('-c','--config',default='irc_eda/config/project.yaml')
    args=parser.parse_args(argv)
    try:
        p=Pipeline(load_config(args.config))
        if args.command=='inspect': result={'points':len(p.inspect()),'inventory_files':len(p.inventory)}
        elif args.command=='prepare':result={'prepared':len(p.prepare())}
        elif args.command in ('run-iqa','run-eda'):result=p.run(args.command.split('-')[1])
        elif args.command=='analyze':result=p.analyze()
        elif args.command=='check-samples':result=samples(p)
        else:
            p.inspect();p.prepare();samples(p);p.run('iqa');p.run('eda');result=p.analyze()
        print(json.dumps(result,ensure_ascii=False,indent=2))
        return 0
    except (ValueError,FileNotFoundError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr);return 2


if __name__=='__main__':raise SystemExit(main())
