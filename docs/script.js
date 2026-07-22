document.addEventListener('DOMContentLoaded', function () {
  const sections = document.querySelectorAll('.variable-section');
  sections.forEach(function (section) {
    const tabs = section.querySelectorAll('.tab');
    const panels = section.querySelectorAll('.panel');
    tabs.forEach(function (tab) {
      tab.addEventListener('click', function () {
        const target = tab.getAttribute('data-projection');
        // Update this section's tabs only — sections are independent,
        // so switching projection on temperature doesn't touch wind, etc.
        tabs.forEach(function (t) {
          const isMatch = t === tab;
          t.classList.toggle('is-active', isMatch);
          t.setAttribute('aria-selected', isMatch ? 'true' : 'false');
        });
        panels.forEach(function (panel) {
          const isMatch = panel.getAttribute('data-projection') === target;
          panel.classList.toggle('is-active', isMatch);
          // Pause any playing video in panels that just became hidden,
          // so switching tabs doesn't leave it running in the background.
          if (!isMatch) {
            const video = panel.querySelector('video');
            if (video && !video.paused) video.pause();
          }
        });
      });
    });
  });
  // Playback speed controls — each panel with a video gets its own
  // Slow/Normal/Fast buttons, scoped so clicking one only affects
  // the video in that same panel.
  const speedButtons = document.querySelectorAll('.speed-btn');
  speedButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      const panel = btn.closest('.panel');
      const video = panel ? panel.querySelector('video') : null;
      if (!video) return;
      video.playbackRate = parseFloat(btn.getAttribute('data-speed'));
      // Only toggle active state among this panel's own speed buttons.
      const siblingButtons = panel.querySelectorAll('.speed-btn');
      siblingButtons.forEach(function (b) {
        b.classList.toggle('is-active', b === btn);
      });
    });
  });
  // Loop toggle — each panel with a video gets its own Loop button,
  // scoped so toggling one only affects the video in that same panel.
  const loopButtons = document.querySelectorAll('[data-loop-toggle]');
  loopButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      const panel = btn.closest('.panel');
      const video = panel ? panel.querySelector('video') : null;
      if (!video) return;
      video.loop = !video.loop;
      btn.classList.toggle('is-active', video.loop);
      btn.setAttribute('aria-pressed', video.loop ? 'true' : 'false');
    });
  });
  // "About" text collapses to 2 lines with a Read more / Show less toggle,
  // scoped per-section so each variable's description is independent.
  const aboutWraps = document.querySelectorAll('.meta-about-wrap');
  aboutWraps.forEach(function (wrap) {
    const text = wrap.querySelector('.meta-about');
    const toggle = wrap.querySelector('.meta-about-toggle');
    if (!text || !toggle) return;
    const isOverflowing = text.scrollHeight > text.clientHeight + 1;
    if (!isOverflowing) {
      toggle.style.display = 'none';
      return;
    }
    toggle.addEventListener('click', function () {
      const expanded = text.classList.toggle('is-expanded');
      toggle.textContent = expanded ? 'Show less' : 'Read more';
    });
  });
});
