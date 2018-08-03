import argparse
import logging
import os

import tornado.escape
import tornado.httpserver
import tornado.ioloop
import tornado.web

from recipe_extractor.recipe_extractor import RecipeExtractor


class RecipeExtractorHandler(tornado.web.RequestHandler):
    route = r'/recipe_extractor'
    version = '2018080200'
    _re = None

    def initialize(self, *args, **kwargs):
        if RecipeExtractorHandler._re is None:
            RecipeExtractorHandler._re = RecipeExtractor()
            self._re = RecipeExtractorHandler._re

    def post(self):
        def error_wrong_format():
            self.set_status(400)
            self.write({
                'status': False,
                'message': 'Input must be list of dictionaries containing "syn_paragraph", "doi", "abstract" fields.'
            })

        documents = tornado.escape.json_decode(self.request.body)

        if not isinstance(documents, list):
            return error_wrong_format()
        for x in documents:
            if not isinstance(x, dict) \
                    or 'syn_paragraph' not in x \
                    or 'doi' not in x \
                    or 'abstract' not in x:
                return error_wrong_format()

        results = []

        for doc in documents:
            structures_list, fails = self._re.get_materials(
                doc['doi'], doc['abstract'], doc['syn_paragraph']
            )
            result = {
                'structures': structures_list,
                'failed': fails,
            }

            results.append(result)

        self.write({
            'status': True,
            'results': results,
            'version': self.version
        })


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.DEBUG)

    parser = argparse.ArgumentParser(description='Synthesis Project Web Service.')
    parser.add_argument('--address', action='store', type=str, default='127.0.0.1',
                        help='Listen address.')
    parser.add_argument('--port', action='store', type=int, default=7731,
                        help='Listen port.')

    args = parser.parse_args()

    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'

    server = tornado.httpserver.HTTPServer(
        tornado.web.Application(
            [
                (RecipeExtractorHandler.route, RecipeExtractorHandler)
            ]
        )
    )
    server.bind(address=args.address, port=args.port)
    logging.info('Going to main loop on %s:%d...', args.address, args.port)
    logging.info('Spawning processes...')
    server.start()
    try:
        tornado.ioloop.IOLoop.current().start()
    except KeyboardInterrupt:
        pass
