import test from 'node:test'
import assert from 'node:assert/strict'
import { clampFigureZoom, figureZoomScroll, fitFigure } from '../src/utils/figureViewport.js'

test('fit includes every edge of wide and tall scientific figures without changing aspect ratio', () => {
  assert.deepEqual(fitFigure(4000, 1000, 800, 600), { width: 800, height: 200 })
  assert.deepEqual(fitFigure(1000, 4000, 800, 600), { width: 150, height: 600 })
  assert.deepEqual(fitFigure(200, 100, 800, 600), { width: 200, height: 100 })
})

test('unloaded images and unavailable viewport dimensions cannot produce invalid CSS sizes', () => {
  for (const dimensions of [[0, 400, 800, 600], [400, 400, 800, 0], [NaN, 400, 800, 600], [400, Infinity, 800, 600]]) {
    assert.deepEqual(fitFigure(...dimensions), { width: 0, height: 0 })
  }
  assert.equal(clampFigureZoom(NaN), 1)
  assert.equal(clampFigureZoom(-2), 1)
  assert.equal(clampFigureZoom(99), 8)
})

test('zooming a centred short figure preserves the image centre and only scrolls overflowing axes', () => {
  assert.deepEqual(figureZoomScroll({ imageWidth: 800, imageHeight: 200, viewportWidth: 800, viewportHeight: 600, nextImageWidth: 1600, nextImageHeight: 400 }), { left: 400, top: 0 })
})

test('zoom preserves the inspected image point after the user has panned away from the centre', () => {
  assert.deepEqual(figureZoomScroll({ scrollLeft: 1000, scrollTop: 300, imageWidth: 2000, imageHeight: 1200, viewportWidth: 800, viewportHeight: 600, nextImageWidth: 3000, nextImageHeight: 1800 }), { left: 1700, top: 600 })
})

test('resetting zoom recentres both axes and cannot leave blank scroll space', () => {
  assert.deepEqual(figureZoomScroll({ scrollLeft: 1200, scrollTop: 600, imageWidth: 2000, imageHeight: 1200, viewportWidth: 800, viewportHeight: 600, nextImageWidth: 800, nextImageHeight: 480 }), { left: 0, top: 0 })
})

test('fullscreen resize preserves the inspected point and clamps image edges to the new viewport', () => {
  assert.deepEqual(figureZoomScroll({ scrollLeft: 1200, scrollTop: 600, imageWidth: 2000, imageHeight: 1200, viewportWidth: 800, viewportHeight: 600, nextImageWidth: 2000, nextImageHeight: 1200, nextViewportWidth: 1200, nextViewportHeight: 900 }), { left: 800, top: 300 })
})
