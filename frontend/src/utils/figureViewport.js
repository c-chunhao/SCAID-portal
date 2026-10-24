export const MIN_FIGURE_ZOOM = 1
export const MAX_FIGURE_ZOOM = 8

const positive = value => Number.isFinite(value) && value > 0
const clamp = (value, minimum, maximum) => Math.min(maximum, Math.max(minimum, value))

// A fit always includes the entire image, without enlarging small originals.
export function fitFigure(width, height, viewportWidth, viewportHeight) {
  if (![width, height, viewportWidth, viewportHeight].every(positive)) return { width: 0, height: 0 }
  const scale = Math.min(1, viewportWidth / width, viewportHeight / height)
  return { width: width * scale, height: height * scale }
}

export function clampFigureZoom(zoom) {
  return Number.isFinite(zoom) ? clamp(zoom, MIN_FIGURE_ZOOM, MAX_FIGURE_ZOOM) : MIN_FIGURE_ZOOM
}

// Keep the same image point under the viewport's centre through zoom or resize.
// Images smaller than the viewport are centred within the scrollable stage.
export function figureZoomScroll({
  scrollLeft = 0, scrollTop = 0,
  imageWidth, imageHeight, viewportWidth, viewportHeight,
  nextImageWidth, nextImageHeight,
  nextViewportWidth = viewportWidth, nextViewportHeight = viewportHeight
}) {
  const axis = (scroll, size, viewport, nextSize, nextViewport) => {
    if (![size, viewport, nextSize, nextViewport].every(positive)) return 0
    const oldInset = Math.max(0, (viewport - size) / 2)
    const newInset = Math.max(0, (nextViewport - nextSize) / 2)
    const centre = clamp((scroll + viewport / 2 - oldInset) / size, 0, 1)
    return clamp(newInset + centre * nextSize - nextViewport / 2, 0, Math.max(0, nextSize - nextViewport))
  }
  return {
    left: axis(scrollLeft, imageWidth, viewportWidth, nextImageWidth, nextViewportWidth),
    top: axis(scrollTop, imageHeight, viewportHeight, nextImageHeight, nextViewportHeight)
  }
}
